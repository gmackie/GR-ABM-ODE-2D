#include "paramwindow.h"
#include "ui_paramwindow.h"
#include "simulation/params.h"
#include "simulation.h"
#include <QStandardItem>
#include <QStandardItemModel>
#include <QMessageBox>
#include <QItemEditorFactory>
#include <QComboBox>
#include <QLineEdit>
#include <QFileDialog>
#include <assert.h>

#include "simulation/xmlhandler.h"
#include "simulation/jsonhandler.h"
#include "simulation/infohandler.h"

/// Private class for adding editable items that directly modify the parameter object
template<typename T>
class EditableTreeItem : public QStandardItem {
  T& _data;
  boost::optional<Range<T> > range;
public:
  EditableTreeItem(T& value, const boost::optional<Range<T> >& _range) : QStandardItem(), _data(value), range(_range) {
    setEditable(true);
    setData(QVariant::fromValue(_data), Qt::EditRole);
  }
  int type() const { return UserType + 1; }
  QStandardItem* clone() const { return new EditableTreeItem<T>(_data, range); }
  void setData ( const QVariant & value, int role = Qt::UserRole + 1 ) {
    if(role == Qt::EditRole || role == Qt::DisplayRole) {
      assert(value.canConvert<T>());
      T tmp = qvariant_cast<T>(value);
      if(!!range && !range->contains(tmp)) { QMessageBox::warning(NULL, "Invalid parameter", "Parameter value out of range", QMessageBox::Ok); return; }
      _data = tmp;
    }
    QStandardItem::setData(value, role);
  }
};

template<>
class EditableTreeItem<std::string> : public QStandardItem {
  std::string& _data;
  boost::optional<Range<std::string> > range;
public:
  EditableTreeItem(std::string& value, const boost::optional<Range<std::string> >& _range) : QStandardItem(), _data(value), range(_range) {
    setEditable(true);
    setData(QVariant::fromValue(QString::fromStdString(_data)), Qt::EditRole);
  }
  int type() const { return UserType + 1; }
  QStandardItem* clone() const { return new EditableTreeItem<std::string>(_data, range); }
  void setData ( const QVariant & value, int role = Qt::UserRole + 1 ) {
    if(role == Qt::EditRole || role == Qt::DisplayRole) {
      _data = value.toString().toStdString();
    }
    QStandardItem::setData(value, role);
  }
};

/// Visitor that visits each parameter property and generates a row with a name, editable value, units, and description
struct TreeVisitor {
  QStandardItem* root;
  TreeVisitor(QStandardItem* _root) : root(_root) { Q_CHECK_PTR(root); }
  /// Find the parent item that matches the path, or create elements of the path until we reach the end of the path
  QStandardItem* findParent(QStandardItem* r, const QStringList::iterator& start, const QStringList::iterator& end) {
    if(start == end) return r;
    for(int i=0;i<r->rowCount();i++) {
      QString childData = r->child(i)->data(Qt::DisplayRole).toString();
      QString top = *start;
      if(childData == top)
        return findParent(r->child(i), start + 1, end);
    }
    r->appendRow(new QStandardItem(*start));
    r->child(r->rowCount()-1)->setEditable(false);
    return findParent(r->child(r->rowCount()-1), start+1, end);
  }

  template<typename T>
  void visit(T& val, Params::ParamDescriptor<T>& desc) {
    QStringList path = QString::fromStdString(desc.path).split('.', QString::SkipEmptyParts);
    QStandardItem* parent = findParent(root, path.begin(), path.end());
    QList<QStandardItem*> items;
    items << new QStandardItem(QString::fromStdString(desc.name));
    items.back()->setEditable(false);
    items << new EditableTreeItem<T>(val, desc.range);
    items << new QStandardItem(QString::fromStdString(desc.units));
    items.back()->setEditable(false);
    items << new QStandardItem(QString::fromStdString(desc.desc));
    items.back()->setEditable(false);
    parent->appendRow(items);
  }
};

ParamWindow::ParamWindow(Simulation& _sim, Params& _params, boost::property_tree::ptree& _pt, QWidget *parent) :
  QWidget(parent),
  params(_params),
  pt(_pt),
  sim(_sim),
  ui(new Ui::ParamWindow)
{
  QItemEditorFactory* factory = new QItemEditorFactory();
  QItemEditorCreatorBase* precisionEditCreator = new QStandardItemEditorCreator<PrecisionEdit>();
  factory->registerEditor(QVariant::Double, precisionEditCreator);
  factory->registerEditor(QVariant::Int, new QStandardItemEditorCreator<IntSpin>());
  factory->registerEditor(QVariant::UInt, new QStandardItemEditorCreator<IntSpin>());
  factory->registerEditor(QVariant::Bool, new QStandardItemEditorCreator<BooleanEditor>());
  factory->registerEditor(QVariant::String, new QStandardItemEditorCreator<QLineEdit>());
  QItemEditorFactory::setDefaultFactory(factory);
  ui->setupUi(this);
  model = new QStandardItemModel(0, 4, this);
  reloadParams(params);
  ui->treeView->setModel(model);
}

void ParamWindow::reloadParams(Params &p) {
  model->clear();
  model->insertColumns(0,4);
  model->setHeaderData(0, Qt::Horizontal, "Name");
  model->setHeaderData(1, Qt::Horizontal, "Value");
  model->setHeaderData(2, Qt::Horizontal, "Units");
  model->setHeaderData(3, Qt::Horizontal, "Description");
  TreeVisitor visitor(model->invisibleRootItem());
  p.visitProperties(visitor);
}

ParamWindow::~ParamWindow()
{
  delete ui;
}

void ParamWindow::on_pushButtonSave_clicked()
{
  QString selected;
  QString fileName = QFileDialog::getSaveFileName(this, tr("Save parameter file"), QString(),
                                                  tr("XML Files (*.xml);;JSON Files (*.json);;INFO Files (*.info)"), &selected);
  ParamFileHandler* handler = NULL;
  if(selected.startsWith("XML")) handler = new XMLHandler("GR");
  else if(selected.startsWith("JSON")) handler = new JSONHandler("GR");
  else if(selected.startsWith("INFO")) handler = new INFOHandler("GR");
  else if(fileName.isNull()) return;  // Cancelled
  else { QMessageBox::critical(this, tr("Error"), tr("Invalid format specified")); return; }
  std::ofstream f(fileName.toStdString().c_str());
  params.save(f, handler, pt);  //No need to lock here, just reading to file
  delete handler;
}

void ParamWindow::on_pushButtonLoad_clicked()
{
  QString selected;
  QString fileName = QFileDialog::getOpenFileName(this, tr("Save parameter file"), QString(),
                                                  tr("XML Files (*.xml);;JSON Files (*.json);;INFO Files (*.info)"), &selected);
  ParamFileHandler* handler = NULL;
  if(selected.startsWith("XML")) handler = new XMLHandler("GR");
  else if(selected.startsWith("JSON")) handler = new JSONHandler("GR");
  else if(selected.startsWith("INFO")) handler = new INFOHandler("GR");
  else if(fileName.isNull()) return;
  else { QMessageBox::critical(this, tr("Error"), tr("Invalid format specified")); return; }
  std::ifstream f(fileName.toStdString().c_str());
  try {
    Params tmp;
    tmp.load(f, handler, pt);
    if(!handler->good())
      throw std::runtime_error("See log for details");
    sim.lock(); //Writer lock!
    sim.modelLock();
    params = tmp;
    sim.unlock();
    sim.modelUnlock();
  } catch(std::exception& e) {
    QMessageBox::critical(this, "Error", QString("Unable to load file:\n%1").arg(e.what()));
  }
  f.close();
  delete handler;
  reloadParams(params);
}
