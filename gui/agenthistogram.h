#ifndef AGENTHISTOGRAM_H
#define AGENTHISTOGRAM_H

#include <QMainWindow>
#include <QVector>
#include "agentswidget.h"
#include "simulation/macrophage.h"

namespace Ui
{
class AgentHistogram;
}

class Simulation;
class QwtPlotCurve;
class QwtPlotItem;
class QProgressBar;

struct MinMaxAgentPropVisitor;
struct HistogramAgentPropVisitor;

class AgentHistogram : public QMainWindow
{
  Q_OBJECT
  friend struct MinMaxAgentPropVisitor;
  friend struct HistogramAgentPropVisitor;

  QProgressBar* pbar;

  size_t nbuckets;
  bool _macFilter[Mac::NSTATES][AgentsWidget::NUM];
  bool _drawTgam;
  bool _drawTcyt;
  bool _drawTreg;
  QVector<QwtPlotCurve*> curves;
  QVector<QVector<double> > xdata;
  QVector<QVector<double> > ydata;
  QVector<double> min;
  QVector<double> max;
  const Simulation& sim;

public:
  explicit AgentHistogram(const Simulation& _sim, QWidget *parent=0);
  ~AgentHistogram();
  /*virtual*/
  void closeEvent(QCloseEvent * e)
  {
    emit closing(false);
    QMainWindow::closeEvent(e);
  }
public slots:
  void updatePlot(int);
  void updatePlot();
  void setAgentFilter(int agentType, int state, int secondaryState, bool enabled);
private slots:
  void on_actionSave_triggered();
  void showCurve(QwtPlotItem*);

signals:
  void closing(bool);

private:
  Ui::AgentHistogram *ui;
};

#endif // AGENTHISTOGRAM_H
