#ifndef GRAPHCONTROLLER_H
#define GRAPHCONTROLLER_H

#include "maininterface.h"
#include <QVector>
#include "graphviewer.h"

struct StatCounterVisitor;
struct StatUpdateVisitor;

class GraphController : public QObject {
  Q_OBJECT

  friend struct StatCounterVisitor;
  friend struct StatUpdateVisitor;
  friend class GraphViewer;

  const size_t _maxSamples;

  QVector<double> _xdata;
  QVector<QVector<double> > _ydata;
  QVector<QString> _names;
  QList<GraphViewer*> _windows;

public:
  
  GraphController(const Stats& stats, size_t maxSamples);
  ~GraphController();

  const QVector<double>& getYdata(size_t i) { return _ydata[i]; }
  const QVector<double>& getXdata() { return _xdata; }
  QString& getName(size_t i) { return _names[i]; }

  void update(const Stats& stats, double t);

public slots:
  void showNewTimeGraph();
  void closeGraph(QWidget* w);

signals:
  void updated();
};

#endif //GRAPHCONTROLLER_H
