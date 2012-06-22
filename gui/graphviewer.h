#ifndef GRAPHVIEWER_H
#define GRAPHVIEWER_H

#include <QMainWindow>
#include <qwt_plot.h>
#include <qwt_plot_curve.h>

namespace Ui {
class GraphViewer;
}

class GraphController;

class GraphViewer : public QMainWindow
{
    Q_OBJECT
    Ui::GraphViewer* ui;
    QVector<QwtPlotCurve*> curves;
    GraphController& ctrlr;
public:
    explicit GraphViewer(GraphController&, QWidget* parent=0);
    ~GraphViewer();
signals:
  void closed(QWidget*);
public slots:
  void updatePlot();
  void showCurve(QwtPlotItem*, bool);
  void closeEvent(QCloseEvent*) { emit closed(this); }
private slots:
  void on_actionSave_triggered();

  void on_actionShow_Legend_toggled(bool arg1);
};

#endif // GRAPHVIEWER_H
