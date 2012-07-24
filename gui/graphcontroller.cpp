#include "gui/graphcontroller.h"
#include "simulation/stat.h"

struct StatCounterVisitor {
  GraphController& _gc;
  StatCounterVisitor(GraphController& gc) : _gc(gc) {}
  template<typename T>
  void visit(const char* name, const T&, const char*) {
    _gc._ydata.push_back(QVector<double>());
    _gc._ydata.back().reserve(_gc._maxSamples);
    _gc._names.push_back(QString(name));
  }
  template<typename DataType, size_t N, typename Enum>
  void visit(const char* name, const boost::array<DataType, N>&, const char*, Enum) {
    for(size_t i=0;i<N;i++) {
      _gc._ydata.push_back(QVector<double>());
      _gc._ydata.back().reserve(_gc._maxSamples);
      std::stringstream ss;
      ss<<name<<'.'<<Enum(i);
      _gc._names.push_back(QString::fromStdString(ss.str()));
    }
  }
};

struct StatUpdateVisitor {
  QVector<QVector<double> >& _ydata;
  size_t iter;
  StatUpdateVisitor(QVector<QVector<double> >& ydata) : _ydata(ydata), iter(0) {}
  template<typename T>
  void visit(const char*, const T& val, const char*) {
    _ydata[iter++].back() = double(val);
  }
  template<typename DataType, size_t N, typename Enum>
  void visit(const char*, const boost::array<DataType, N>& val, const char*, Enum) {
    for(size_t i=0;i<N;i++)
      _ydata[iter++].back() = double(val[i]);
  }
};

GraphController::GraphController(const Stats& stats, size_t maxSamples) : _maxSamples(maxSamples) {
   _xdata.reserve(maxSamples);
  StatCounterVisitor visitor(*this);
  stats.visit(visitor);
}

void GraphController::update(const Stats& stats, double t) {
  t /= 144.0; //Convert to days
  if((size_t)_xdata.size() >= _maxSamples)
  {
    memmove(_xdata.data(), _xdata.data()+1, (_xdata.size() - 1)*sizeof(_xdata[0]));
    for(int i=0;i<_ydata.size(); i++)
      memmove(_ydata[i].data(), _ydata[i].data()+1, (_ydata[i].size() - 1)*sizeof(_ydata[i][0]));
  }
  else {
    _xdata.push_back(t);
    for(int i=0;i<_ydata.size(); i++)
      _ydata[i].push_back(0);
  }
  _xdata.back() = t;
  StatUpdateVisitor visitor(_ydata);
  stats.visit(visitor);
  emit updated();
}

void GraphController::showNewTimeGraph() {
  GraphViewer* graphViewer = new GraphViewer(*this);
  connect(this, SIGNAL(updated(void)), graphViewer, SLOT(updatePlot(void)));
  connect(graphViewer, SIGNAL(closed(QWidget*)), this, SLOT(closeGraph(QWidget*)));
  if(!_xdata.isEmpty())
    graphViewer->updatePlot();
  graphViewer->show();
  _windows.push_back(graphViewer);
}

GraphController::~GraphController() {
  for(int i=0;i<_windows.size();i++)
    _windows[i]->deleteLater();
}

void GraphController::closeGraph(QWidget* w) {
  for(int i=0;i<_windows.size();i++)
    if(_windows[i] == w)
    {
      _windows.removeAt(i);
      w->deleteLater();
      break;
    }
}
