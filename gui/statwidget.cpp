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
  StatWidgetNameVisitor(StatWidget& sw) : _sw(sw), iter(0) {}
  template<typename T>
  void visit(const char* name, const T& val, const char* desc)
  {
    _sw._ui.tableWidget->insertRow(iter);
    QTableWidgetItem* item = new QTableWidgetItem();
    item->setText(QString(name));
    item->setToolTip(QString(desc));
    _sw._ui.tableWidget->setItem(iter, 0, item);
    _sw._ui.tableWidget->setItem(iter, 1, new QTableWidgetItem(toString(val)));
    ++iter;
  }
  template<typename T, size_t N, typename E>
  void visit(const char* name, const boost::array<T, N>& val, const char* desc, E e)
  {
    _sw._ui.tableWidget->insertRow(iter);
    QTableWidgetItem* item = new QTableWidgetItem();
    item->setText(QString(name));
    item->setToolTip(QString(desc));
    _sw._ui.tableWidget->setItem(iter, 0, item);
    item = new QTableWidgetItem();
    item->setText(toString(val));
    item->setToolTip(enumstoString(e));
    _sw._ui.tableWidget->setItem(iter, 1, item);
    ++iter;
  }
};

struct StatWidgetUpdateVisitor
{
  StatWidget& _sw;
  int iter;
  StatWidgetUpdateVisitor(StatWidget& sw) : _sw(sw), iter(0) {}
  template<typename T>
  void visit(const char*, const T& val, const char*)
  {
    QTableWidgetItem* item = _sw._ui.tableWidget->item(iter, 1);
    Q_CHECK_PTR(item);
    item->setText(toString(val));
    ++iter;
  }
  template<typename T, size_t N, typename E>
  void visit(const char*, const boost::array<T, N>& val, const char*, E)
  {
    QTableWidgetItem* item = _sw._ui.tableWidget->item(iter, 1);
    Q_CHECK_PTR(item);
    item->setText(toString(val));
    ++iter;
  }
};

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
