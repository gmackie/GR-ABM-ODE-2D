#include "agenthistogram.h"
#include "ui_agenthistogram.h"
#include "simulation.h"
#include <QSpinBox>
#include <QProgressBar>

#include <qwt_plot_curve.h>

#include <QFileDialog>
#include <qwt_legend.h>
#include <qwt_legend_item.h>
#if (QWT_VERSION >> 16) == 0x06
#include <qwt_plot_picker.h>
#include <qwt_plot_renderer.h>
#include <qwt_symbol.h>
#endif

#include <limits>

AgentHistogram::AgentHistogram(const Simulation& _sim, QWidget *parent) :
  QMainWindow(parent),
  nbuckets(10),
  _drawTcyt(true),
  _drawTgam(true),
  _drawTreg(true),
  sim(_sim),
  ui(new Ui::AgentHistogram)
{
  memset(_macFilter, 1, sizeof(_macFilter));
  ui->setupUi(this);
  {
    //Qtcreator doesn't allow other widgets in toolbar, have to add them manually :(
    QSpinBox* spinBucket = new QSpinBox();
    spinBucket->setSingleStep(5);
    spinBucket->setRange(1, std::numeric_limits<int>::max());
    spinBucket->setSuffix(" buckets");
    spinBucket->setValue(nbuckets);
    connect(spinBucket, SIGNAL(valueChanged(int)), SLOT(updatePlot(int)));
    ui->toolBar->addWidget(spinBucket)->setVisible(true);
  }
  {
    pbar = new QProgressBar(this);
    pbar->setValue(100);
    ui->toolBar->addWidget(pbar)->setVisible(true);
  }
#if (QWT_VERSION >> 16) == 0x06
  new QwtPlotPicker(QwtPlot::xBottom, QwtPlot::yLeft, QwtPicker::CrossRubberBand,
                    QwtPlotPicker::AlwaysOn, ui->centralwidget->canvas());
#endif

  QwtLegend* _legend = new QwtLegend();
  _legend->setItemMode(QwtLegend::ClickableItem);
  ui->centralwidget->insertLegend(_legend, QwtPlot::RightLegend);
  connect(ui->centralwidget, SIGNAL(legendClicked(QwtPlotItem*)), SLOT(showCurve(QwtPlotItem*)));
}

AgentHistogram::~AgentHistogram()
{
  delete ui;
}

void AgentHistogram::updatePlot(int sz)
{
  nbuckets = sz;
  updatePlot();
}

struct MinMaxAgentPropVisitor
{
  AgentHistogram& ah;
  int iter;
  MinMaxAgentPropVisitor(AgentHistogram& _ah) : ah(_ah), iter(0) {}
  void visit(const Agent* a)
  {
    switch(a->getAgentType())
      {
      case MAC:
      {
        const Mac* pMac = static_cast<const Mac*>(a);
        Q_CHECK_PTR(pMac);
        char state = (pMac->getNFkB() << AgentsWidget::NFKB)
                     | (pMac->getStat1() << AgentsWidget::STAT1)
                     | (pMac->isDeactivated() << AgentsWidget::DEACT);
        state |= (state == 0) << AgentsWidget::OTHER;
        state &= (ah._macFilter[pMac->getState()][AgentsWidget::NFKB] << AgentsWidget::NFKB)
                 | (ah._macFilter[pMac->getState()][AgentsWidget::STAT1] << AgentsWidget::STAT1)
                 | (ah._macFilter[pMac->getState()][AgentsWidget::DEACT] << AgentsWidget::DEACT)
                 | (ah._macFilter[pMac->getState()][AgentsWidget::OTHER] << AgentsWidget::OTHER);
        if(state && ah._macFilter[pMac->getState()][AgentsWidget::ENBL])
          a->visitProperties(*this);
        break;
      }
      case TGAM:
        if(ah._drawTgam) a->visitProperties(*this);
        break;
      case TREG:
        if(ah._drawTgam) a->visitProperties(*this);
        break;
      case TCYT:
        if(ah._drawTgam) a->visitProperties(*this);
        break;
      default:
        break;
      }
  }
  template<typename T> void visit(const char*, const T&, const char*) {}  //Skip the properties we can't deal with
};

#ifndef unlikely //branch optimization
#define unlikely(x) (x)
#endif

template<>
void MinMaxAgentPropVisitor::visit(const char *name, const Scalar &val, const char *desc)
{
  if(unlikely(ah.curves.size() <= iter))
    {
      ah.curves.append(new QwtPlotCurve(QString(name)));
      ah.curves[iter]->setStyle(QwtPlotCurve::Sticks);    //Sticks? or Steps?
      ah.curves[iter]->attach(ah.ui->centralwidget);
      ah.curves[iter]->setVisible(false);
      ah.curves[iter]->legendItem()->setToolTip(desc);
      ah.min.append(std::numeric_limits<Scalar>::max());  //These are reversed in order to force taking a real value
      ah.max.append(-std::numeric_limits<Scalar>::max());
      ah.xdata.append(QVector<double>(ah.nbuckets));
      ah.ydata.append(QVector<double>(ah.nbuckets));
    }
  ah.min[iter] = std::min(ah.min[iter], val);
  ah.max[iter] = std::max(ah.max[iter], val);
  ++iter;
}

struct HistogramAgentPropVisitor
{
  AgentHistogram& ah;
  int iter;
  HistogramAgentPropVisitor(AgentHistogram& _ah) : ah(_ah), iter(0) {}
  void visit(const Agent* a)
  {
    switch(a->getAgentType())
      {
      case MAC:
      {
        const Mac* pMac = static_cast<const Mac*>(a);
        Q_CHECK_PTR(pMac);
        char state = (pMac->getNFkB() << AgentsWidget::NFKB)
                     | (pMac->getStat1() << AgentsWidget::STAT1)
                     | (pMac->isDeactivated() << AgentsWidget::DEACT);
        state |= (state == 0) << AgentsWidget::OTHER;
        state &= (ah._macFilter[pMac->getState()][AgentsWidget::NFKB] << AgentsWidget::NFKB)
                 | (ah._macFilter[pMac->getState()][AgentsWidget::STAT1] << AgentsWidget::STAT1)
                 | (ah._macFilter[pMac->getState()][AgentsWidget::DEACT] << AgentsWidget::DEACT)
                 | (ah._macFilter[pMac->getState()][AgentsWidget::OTHER] << AgentsWidget::OTHER);
        if(state && ah._macFilter[pMac->getState()][AgentsWidget::ENBL])
          a->visitProperties(*this);
        break;
      }
      case TGAM:
        if(ah._drawTgam) a->visitProperties(*this);
        break;
      case TREG:
        if(ah._drawTgam) a->visitProperties(*this);
        break;
      case TCYT:
        if(ah._drawTgam) a->visitProperties(*this);
        break;
      default:
        break;
      }
  }
  template<typename T> void visit(const char*, const T&, const char*) {}  //Skip the properties we can't deal with
};

template<>
void HistogramAgentPropVisitor::visit(const char *, const Scalar & val, const char *)
{
  int index = int((val-ah.min[iter]) * (ah.nbuckets-1) / (ah.max[iter] - ah.min[iter]));
  ++ah.ydata[iter][index];
  ++iter;
}

void AgentHistogram::updatePlot()
{
  pbar->setValue(0);
  pbar->repaint();
  sim.lock(); //Locking here because this could be called without a lock (via spinBucket)
  {
    this->min.fill(std::numeric_limits<double>::max());
    this->max.fill(-std::numeric_limits<double>::max());
    MinMaxAgentPropVisitor visitor(*this);
    for(size_t i=0; i<sim.getMacList().size(); i++)
      {
        visitor.iter = 0;
        visitor.visit(&(sim.getMacList()[i]));
      }
    pbar->setValue(10);
    pbar->repaint();
    for(size_t i=0; i<sim.getTgamList().size(); i++)
      {
        visitor.iter = 0;
        visitor.visit(&(sim.getTgamList()[i]));
      }
    pbar->setValue(20);
    pbar->repaint();
    for(size_t i=0; i<sim.getTregList().size(); i++)
      {
        visitor.iter = 0;
        visitor.visit(&(sim.getTregList()[i]));
      }
    pbar->setValue(30);
    pbar->repaint();
    for(size_t i=0; i<sim.getTcytList().size(); i++)
      {
        visitor.iter = 0;
        visitor.visit(&(sim.getTcytList()[i]));
      }
    pbar->setValue(40);
    pbar->repaint();
  }
  for(int i=0; i<curves.size(); i++)
    {
      if(unlikely(xdata[i].size() != (int)nbuckets))
        {
          xdata[i].resize(nbuckets);    //resize() will do nothing if vector is the correct size
          ydata[i].resize(nbuckets);
        }
#if (QWT_VERSION >> 16) == 5
      curves[i]->setRawData(xdata[i].data(), ydata[i].data(), nbuckets);
#else
      curves[i]->setRawSamples(xdata[i].data(), ydata[i].data(), nbuckets);
#endif
      ydata[i].fill(0);
      if(this->min[i] >= this->max[i])     //To prevent divide by zero, double the width
        {
          this->max[i] =  this->min[i]+1;
        }
      for(size_t j=0; j<nbuckets; j++)
        xdata[i][j] = (j * (this->max[i] - this->min[i]) / nbuckets) + this->min[i];
    }
  pbar->setValue(50);
  pbar->repaint();
  {
    HistogramAgentPropVisitor visitor(*this);
    for(size_t i=0; i<sim.getMacList().size(); i++)
      {
        visitor.iter = 0;
        visitor.visit(&(sim.getMacList()[i]));
      }
    pbar->setValue(60);
    pbar->repaint();
    for(size_t i=0; i<sim.getTgamList().size(); i++)
      {
        visitor.iter = 0;
        visitor.visit(&(sim.getTgamList()[i]));
      }
    pbar->setValue(70);
    pbar->repaint();
    for(size_t i=0; i<sim.getTregList().size(); i++)
      {
        visitor.iter = 0;
        visitor.visit(&(sim.getTregList()[i]));
      }
    pbar->setValue(80);
    pbar->repaint();
    for(size_t i=0; i<sim.getTcytList().size(); i++)
      {
        visitor.iter = 0;
        visitor.visit(&(sim.getTcytList()[i]));
      }
    pbar->setValue(90);
    pbar->repaint();
  }
  sim.unlock();
  ui->centralwidget->replot();
  pbar->setValue(100);
  pbar->repaint();
}

void AgentHistogram::on_actionSave_triggered()
{
  QString fname = QFileDialog::getSaveFileName(this, tr("Save As"), ".", tr("Images (*.png, *.xpm, *.jpg)"));
#if (QWT_VERSION >> 16) == 0x05
  QImage img(width(), height(), QImage::Format_RGB32);
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

void AgentHistogram::showCurve(QwtPlotItem* curve)
{
  curve->show();
  //Disable everyone else
  for(int i=0; i<curves.size(); i++)
    {
      QwtPlotCurve* c = curves[i];
      if(c != curve)
        c->hide();
    }
  ui->centralwidget->setTitle(curve->title());
  ui->centralwidget->replot();
}

void AgentHistogram::setAgentFilter(int agentType, int state, int secondaryState, bool enabled)
{
  switch((AgentType)agentType)
    {
    case MAC:
      _macFilter[state][secondaryState] = enabled;
      break;
    case TGAM:
      _drawTgam = secondaryState == AgentsWidget::ENBL && enabled;
      break;
    case TCYT:
      _drawTcyt = secondaryState == AgentsWidget::ENBL && enabled;
      break;
    case TREG:
      _drawTreg = secondaryState == AgentsWidget::ENBL && enabled;
      break;
    default:
      break;
    }
  if(!isHidden())
    updatePlot();
}
