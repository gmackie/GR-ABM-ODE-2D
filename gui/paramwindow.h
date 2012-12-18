#ifndef PARAMWINDOW_H
#define PARAMWINDOW_H

#include <QSpinBox>
#include <QLineEdit>
#include <QComboBox>
#include <limits>
#include <boost/property_tree/ptree.hpp>

namespace Ui {
  class ParamWindow;
}

class Params;
class Simulation;
class QStandardItemModel;

class IntSpin : public QSpinBox {
  Q_OBJECT
public:
  IntSpin(QWidget* parent=0) : QSpinBox(parent) {
    setMinimum(std::numeric_limits<int>::min());
    setMaximum(std::numeric_limits<int>::max());
  }
};

class PrecisionEdit : public QLineEdit {
  Q_OBJECT
  Q_PROPERTY(double value READ value WRITE setValue USER true)
public:
  PrecisionEdit(QWidget* parent=0) : QLineEdit(parent) {
    this->setValidator(new QDoubleValidator());
  }
  virtual ~PrecisionEdit() {}
  double value() const { return text().toDouble(); }
  void setValue(double v) { setText(QLocale::system().toString(v, 'g', 100)); }
};

class BooleanEditor : public QComboBox {
  Q_OBJECT
  Q_PROPERTY(bool value READ value WRITE setValue USER true)
public:
  BooleanEditor(QWidget* parent=0) : QComboBox(parent) {
    this->addItems(QStringList("false") << "true");

  }
  virtual ~BooleanEditor() {}
  bool value() const { return (bool)this->currentIndex(); }
  void setValue(bool v) { setCurrentIndex(v); }
};

class ParamWindow : public QWidget
{
  Q_OBJECT
  Params& params;
  boost::property_tree::ptree& pt;
  Simulation& sim;
  QStandardItemModel* model;
  
public:
  ParamWindow(Simulation& _sim, Params& _params, boost::property_tree::ptree& _pt, QWidget *parent = 0);
  ~ParamWindow();
  void reloadParams(Params& p);
  
private slots:
  void on_pushButtonSave_clicked();

  void on_pushButtonLoad_clicked();

private:
  Ui::ParamWindow *ui;
};

#endif // PARAMWINDOW_H
