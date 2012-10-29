#include "statwidget.h"
#include "simulation/stat.h"


template<typename T>
QString toString(const T& v)
{
  return QString("%1").arg(v);
}
template<>
QString toString<Scalar>(const Scalar& v)
{
  return QString::number(v, 'f', 2);
}

template<typename T, size_t N>
static QString toString(const boost::array<T, N>& arr)
{
  QString str('(');
  for(size_t i=0; i<N-1; i++)
    str += toString(arr[i]) + ',';
  str += toString(arr[N-1]) + ')';
  return str;
}

template<typename E>
static QString enumstoString(E final)
{
  std::stringstream ss;
  int i=0;
  for(i=0; i<final-1; i++)
    ss<<E(i)<<',';
  ss<<E(i);
  return QString::fromStdString(ss.str());
}

struct StatWidgetNameVisitor
{
  StatWidget& _sw;
  int iter;
  QTreeWidgetItem* root;
  StatWidgetNameVisitor(StatWidget& sw) : _sw(sw), iter(0), root(NULL) {}
  template<typename T>
  void visit(const char* name, const T& val, const char* desc)
  {
    QTreeWidgetItem* item = new QTreeWidgetItem();
    item->setText(0, QString(name));
    item->setText(1, toString(val));
    item->setToolTip(0, QString(desc));
    if(root == NULL)
        _sw._ui.treeWidget->addTopLevelItem(item);
    else
        root->addChild(item);
  }
  template<typename T, size_t N, typename E>
  void visit(const char* name, const boost::array<T, N>& val, const char* desc, E e)
  {
    QTreeWidgetItem* item = new QTreeWidgetItem();
    item->setText(0, QString(name));
    item->setText(1, toString(val));
    item->setToolTip(0, toString(desc));
    item->setToolTip(1, enumstoString(e));
    if(root == NULL)
      _sw._ui.treeWidget->addTopLevelItem(item);
    else
      root->addChild(item);
  }
};

template<>
void StatWidgetNameVisitor::visit(const char* name, const Stats::group_type&, const char* desc) {
    root = new QTreeWidgetItem();
    root->setText(0, QString(name));
    root->setText(1, QString(desc));
    _sw._ui.treeWidget->addTopLevelItem(root);
}

struct StatWidgetUpdateVisitor
{
  StatWidget& _sw;
  int iter;
  QTreeWidgetItem* root;
  StatWidgetUpdateVisitor(StatWidget& sw) : _sw(sw), iter(0), root(NULL) {}
  template<typename T>
  void visit(const char* name, const T& val, const char*)
  {
    QTreeWidgetItem* item = NULL;
    if(!root)
      item = _sw._ui.treeWidget->findItems(name, Qt::MatchExactly, 0).first();
    else
      item = root->child(iter++);
    Q_CHECK_PTR(item);
    item->setText(1, toString(val));
  }
  template<typename T, size_t N, typename E>
  void visit(const char* name, const boost::array<T, N>& val, const char*, E)
  {
    QTreeWidgetItem* item = NULL;
    if(!root)
      item = _sw._ui.treeWidget->findItems(name, Qt::MatchExactly, 0).first();
    else
      item = root->child(iter++);
    Q_CHECK_PTR(item);
    item->setText(1, toString(val));
  }
};

template<>
void StatWidgetUpdateVisitor::visit(const char * name, const Stats::group_type &, const char *)
{   //Sets up the next group to look in
    QString n(name);
    size_t sz = _sw._ui.treeWidget->topLevelItemCount();
    for(size_t i=0;i<sz;i++)
      if(n==_sw._ui.treeWidget->topLevelItem(i)->text(0))
        root =_sw._ui.treeWidget->topLevelItem(i);
    iter = 0;
}

StatWidget::StatWidget(const Stats& stats, QWidget* pParent)
  : QWidget(pParent)
{
  _ui.setupUi(this);
  StatWidgetNameVisitor visitor(*this);
  stats.visit(visitor);
  update();
}

StatWidget::~StatWidget()
{
}

void StatWidget::updateLabels(const Stats& stats)
{
  StatWidgetUpdateVisitor visitor(*this);
  stats.visit(visitor);
}
