#include "gui/graphviewer.h"
#include "gui/graphcontroller.h"
#include "ui_graphviewer.h"

#include <QFileDialog>
#include <qwt_legend.h>
#include <qwt_legend_item.h>
#if (QWT_VERSION >> 16) == 0x06
	#include <qwt_plot_renderer.h>
  #include <qwt_symbol.h>
#endif

GraphViewer::GraphViewer(GraphController& _ctrlr, QWidget* parent) : QMainWindow(parent), ui(new Ui::GraphViewer), ctrlr(_ctrlr) {
  ui->setupUi(this);

  ui->centralwidget->setGeometry(0,0,640,480);

  QwtLegend* _legend = new QwtLegend();
  _legend->setItemMode(QwtLegend::CheckableItem);
  ui->centralwidget->insertLegend(_legend, QwtPlot::RightLegend);

  ui->centralwidget->setAxisTitle(QwtPlot::xBottom, "Time, t, in days");

  static const Qt::GlobalColor validColors[] = {
      Qt::black,
      Qt::red,
      Qt::darkRed,
      Qt::green,
      Qt::darkGreen,
      Qt::blue,
      Qt::darkBlue,
      Qt::cyan,
      Qt::darkCyan,
      Qt::magenta,
      Qt::darkMagenta,
      Qt::darkYellow
  };
  static const size_t ncolors = sizeof(validColors) / sizeof(Qt::GlobalColor);

  curves.resize(ctrlr._names.size());
  for(int i=0;i<ctrlr._names.size();i++) {
    curves[i] = new QwtPlotCurve();
    curves[i]->attach(ui->centralwidget);
    curves[i]->setPen(QPen(QColor(validColors[i%ncolors])));
#if (QWT_VERSION >> 16) == 0x06
    QwtSymbol* sym = new QwtSymbol();
    sym->setStyle(QwtSymbol::Style((i/ncolors) % ((int)QwtSymbol::Hexagon + 1) - 1));
    sym->setSize(QSize(2, 2));
    sym->setColor(QColor(validColors[i%ncolors]));
    sym->setPen(QPen(QColor(validColors[i%ncolors])));
    curves[i]->setSymbol(sym);
    curves[i]->setLegendAttribute(QwtPlotCurve::LegendShowLine);
    curves[i]->setLegendAttribute(QwtPlotCurve::LegendShowSymbol);
#endif
    curves[i]->setTitle(ctrlr._names[i]);
    this->showCurve(curves[i], false);
  }
  connect(ui->centralwidget, SIGNAL(legendChecked(QwtPlotItem*,bool)), SLOT(showCurve(QwtPlotItem*,bool)));
}

void GraphViewer::updatePlot() {
  ui->centralwidget->setAxisScale(QwtPlot::xBottom, ctrlr.getXdata().front(), ctrlr.getXdata().back());
  for(int i=0;i<curves.size();i++)
#if (QWT_VERSION >> 16) == 5
    curves[i]->setRawData(ctrlr.getXdata().data(), ctrlr.getYdata(i).data(), ctrlr.getXdata().size());
#else
    curves[i]->setRawSamples(ctrlr.getXdata().data(), ctrlr.getYdata(i).data(), ctrlr.getXdata().size());
#endif
  ui->centralwidget->replot();
}

GraphViewer::~GraphViewer() {
  for(int i=0;i<curves.size();i++)
    delete curves[i];
  delete ui;
}

void GraphViewer::showCurve(QwtPlotItem* item, bool enabled) {
  item->setVisible(enabled);
  item->setItemAttribute(QwtPlotItem::AutoScale, enabled);
  QWidget* w = item->legendItem();
  if(w && w->inherits("QwtLegendItem")) {
    static_cast<QwtLegendItem*>(w)->setItemMode(QwtLegend::CheckableItem);
    static_cast<QwtLegendItem*>(w)->setChecked(enabled);
  }
  ui->centralwidget->replot();
}

void GraphViewer::on_actionSave_triggered()
{
    QString fname = QFileDialog::getSaveFileName(this, tr("Save As"), ".", tr("Images (*.png, *.xpm, *.jpg)"));
#if (QWT_VERSION >> 16) == 0x05
    QImage img(640, 480, QImage::Format_RGB32);
    img.fill(Qt::white);
    QwtPlotPrintFilter filter;
    filter.setOptions(QwtPlotPrintFilter::PrintFrameWithScales | QwtPlotPrintFilter::PrintBackground);
    ui->centralwidget->print(img, filter);
    img.save(fname);
#else	//QWT_VERSION == 0x06
    QwtPlotRenderer renderer(this);
    renderer.setDiscardFlags(QwtPlotRenderer::DiscardLegend);
    renderer.renderDocument(ui->centralwidget, fname, ui->centralwidget->size());
#endif
}

void GraphViewer::on_actionShow_Legend_toggled(bool )
{
    //ui->centralwidget->insertLegend(arg1 ? _legend : NULL, QwtPlot::ExternalLegend);
}
