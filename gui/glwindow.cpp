#include "glwindow.h"
#include "colormaps/colormap.h"
#include "maininterface.h"
#include "scalardatasets/scalaragentgrid.h"
#include <QTreeWidget>
#include <QTreeWidgetItemIterator>
#include <QImage>
#include <QCloseEvent>

GLWindow::GLWindow(MainInterface* pItfc, QWidget* parent)
  : QWidget(parent)
  , _ui()
  , _pItfc(pItfc)
  , _selRow(-1)
  , _selCol(-1)
  , _printTime(_PRINT_TIME)
  , _printOutcome(_PRINT_OUTCOME)
  , _trackid(-1)
{
  for (int i = 0; i < NOUTCOMES; i++)
    {
      _status[i] = GR_NONE;
    }

  _ui.setupUi(this);

  _ui.glWidget->dim = pItfc->getSimulation().getSize();
  connect(_ui.glWidget, SIGNAL(updateSelection(int, int)), this, SIGNAL(updateSelection(int, int)));
  connect(_ui.glWidget, SIGNAL(visualize(void)), this, SLOT(visualize(void)));
  connect(_ui.glWidget, SIGNAL(printText(void)), this, SLOT(printText(void)));
  connect(this, SIGNAL(set2DView(void)), _ui.glWidget, SLOT(set2DView(void)));
  connect(this, SIGNAL(set3DView(void)), _ui.glWidget, SLOT(set3DViewHeight(void)));

  QStringList hdrs;
  hdrs << "Property" << "Value" << "Description";
  agentInfoWindow = new QTreeWidget();
  agentInfoWindow->resize(400,400);
  agentInfoWindow->setHeaderLabels(hdrs);
  agentInfoWindow->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
  agentInfoWindow->setWindowTitle("Agent Information");
  connect(agentInfoWindow, SIGNAL(itemChanged(QTreeWidgetItem*,int)), SLOT(updateTracking(QTreeWidgetItem*)));
}

GLWindow::~GLWindow()
{
  agentInfoWindow->close();
  delete agentInfoWindow;
}

void GLWindow::updateTracking(QTreeWidgetItem* item)
{
  if(item->checkState(0) == Qt::Unchecked && _trackid == item->data(1, Qt::DisplayRole).toInt())
    {
      _trackid = -1;
      return;
    } //otherwise...
  _trackid = item->data(1, Qt::DisplayRole).toInt();
  agentInfoWindow->blockSignals(true);    //Don't generate any more signals until we update the rest of the gui

  for(QTreeWidgetItemIterator it(agentInfoWindow); *it; it++)
    {
      if((*it) == item) continue;
      (*it)->setCheckState(0, Qt::Unchecked);
    }
  agentInfoWindow->blockSignals(false);
  updateSelectedCellStats();
}

void GLWindow::setColorMap(ColorMap* pCurrentColorMap)
{
  _ui.colorMapWidget->setColorMap(pCurrentColorMap);
}

void GLWindow::updateWindow()
{
  _ui.glWidget->updateGL();
}

void GLWindow::toggleFullScreen()
{
  if (isFullScreen())
    {
      _ui.colorMapWidget->show();
      _ui.labelLegend1->show();
      _ui.labelLegend2->show();
      _ui.labelLegend3->show();
      _ui.labelLegend4->show();
      _ui.labelLegend5->show();

      showNormal();
    }
  else
    {
      _ui.colorMapWidget->hide();
      _ui.labelLegend1->hide();
      _ui.labelLegend2->hide();
      _ui.labelLegend3->hide();
      _ui.labelLegend4->hide();
      _ui.labelLegend5->hide();

      showFullScreen();
    }
}

template<typename T>
static inline QTreeWidgetItem* make_item(const char* name, const T& val, const char* desc)
{
  QTreeWidgetItem* item = new QTreeWidgetItem();
  item->setText(0, QString(name));
  item->setData(1, Qt::DisplayRole, QVariant::fromValue(val));
  item->setText(2, QString(desc));
  return item;
}

template<>
inline QTreeWidgetItem* make_item(const char* name, const Pos& val, const char* desc)
{
  std::stringstream ss;
  ss<<val;
  QTreeWidgetItem* item = new QTreeWidgetItem();
  item->setText(0, QString(name));
  item->setText(1, QString::fromStdString(ss.str()));
  item->setText(2, QString(desc));
  return item;
}

struct AgentInfoVisitor
{
  QTreeWidget* _view;
  QTreeWidgetItem* _curItem;
  const int _trackid;
  AgentInfoVisitor(QTreeWidget* view, int trackid) : _view(view), _curItem(NULL), _trackid(trackid) {}
  void visit(const Agent* a)
  {
    if(!a) return;
    std::stringstream ss;
    _curItem = new QTreeWidgetItem();
    _view->addTopLevelItem(_curItem);
    switch(a->getAgentType())
      {
      case MAC:
        ss<<"Mac "<<(Mac::State)(a->getState());
        visit(static_cast<const Mac*>(a));
        break;
      case TGAM:
        ss<<"Tgam "<<(Tgam::State)(a->getState());
        visit(static_cast<const Tgam*>(a));
        break;
      case TCYT:
        ss<<"Tcyt "<<(Tcyt::State)(a->getState());
        visit(static_cast<const Tcyt*>(a));
        break;
      case TREG:
        ss<<"Treg "<<(Treg::State)(a->getState());
        visit(static_cast<const Treg*>(a));
        break;
      }
    //a->visitProperties(*this);
    _curItem->setText(0, QString::fromStdString(ss.str()));
    _curItem->setData(1, Qt::DisplayRole, QVariant::fromValue(a->getid()));
    _curItem->setExpanded(true);
    _curItem->setFlags(Qt::ItemIsSelectable | Qt::ItemIsUserCheckable | Qt::ItemIsEnabled);
    _curItem->setCheckState(0, (a->getid() == _trackid ? Qt::Checked : Qt::Unchecked));
  }

  void visit(const Mac* a)
  {
    a->visitProperties(*this);
  }
  void visit(const Tgam* a)
  {
    a->visitProperties(*this);
  }
  void visit(const Tcyt* a)
  {
    a->visitProperties(*this);
  }
  void visit(const Treg* a)
  {
    a->visitProperties(*this);
  }

  template<typename T>
  void visit(const char* name, const T& val, const char* desc)
  {
    Q_CHECK_PTR(_curItem);
    _curItem->addChild(make_item(name, val, desc));
  }
};

template<typename AgentType>
static const Agent* findAgent(const std::vector<AgentType>& v, size_t id)
{
  for(size_t i=0; i<v.size(); i++)
    {
      const Agent* a = &(v[i]);
      if(a->getid() == id)
        return a;
    }
  return NULL;
}

void GLWindow::updateSelectedCellStats()
{
  if(_selRow != -1 && _selCol != -1)
    {
      _pItfc->getSimulation().lock();
      if(_trackid > -1)   //Follow the id...
        {
          const Agent* a = findAgent(_pItfc->getSimulation().getMacList(), _trackid);
          if(!a)
            a = findAgent(_pItfc->getSimulation().getTgamList(), _trackid);
          if(!a)
            a = findAgent(_pItfc->getSimulation().getTcytList(), _trackid);
          if(!a)
            a = findAgent(_pItfc->getSimulation().getTregList(), _trackid);
          if(a && (a->getPosition().x != _selRow ||
                   a->getPosition().y != _selCol))
            {
              //Moved, update everything
              emit updateSelection(a->getPosition().x, a->getPosition().y);
            }
        }


      bool topLevelExpanded[3];
      for(size_t i=0; i<3; i++) //Hacky way to keep expanded information
        topLevelExpanded[i] = (agentInfoWindow->topLevelItem(i) ? agentInfoWindow->topLevelItem(i)->isExpanded() : false);
      agentInfoWindow->blockSignals(true);    //Don't generate any more signals until we update the rest of the gui
      agentInfoWindow->clear();   //Not the most efficient
      const ScalarAgentGrid* pAgentGrid = static_cast<ScalarAgentGrid*>(_pItfc->getScalarAgentGrid());
      const ScalarAgentItem& item
        = pAgentGrid->getGrid()[_selRow * _ui.glWidget->dim.x + _selCol];
      QTreeWidgetItem* gridInfo = new QTreeWidgetItem();
      agentInfoWindow->addTopLevelItem(gridInfo);
      gridInfo->setText(0, "Grid Properties");
      gridInfo->addChild(make_item("nKillings",     item._nKillings,     ""));
      gridInfo->addChild(make_item("nRecruitments", item._nRecruitments, ""));
      gridInfo->addChild(make_item("nRecruitmentsMac", item._nRecruitmentsMac, ""));
      gridInfo->addChild(make_item("nRecruitmentsTgam", item._nRecruitmentsTgam, ""));
      gridInfo->addChild(make_item("nRecruitmentsTcyt", item._nRecruitmentsTcyt, ""));
      gridInfo->addChild(make_item("nRecruitmentsTreg", item._nRecruitmentsTreg, ""));
      gridInfo->addChild(make_item("nSecretions",   item._nSecretions,   ""));
      gridInfo->addChild(make_item("attractant",    item._attractant,    ""));
      gridInfo->addChild(make_item("TNF",           item._TNF,           ""));
      gridInfo->addChild(make_item("CCL2",          item._CCL2,          ""));
      gridInfo->addChild(make_item("CCL5",          item._CCL5,          ""));
      gridInfo->addChild(make_item("CXCL9",         item._CXCL9,         ""));
      gridInfo->addChild(make_item("shedTNFR2",     item._shedTNFR2,     ""));
      gridInfo->addChild(make_item("il10",          item._il10,          ""));
      gridInfo->addChild(make_item("extMtb",        item._extMtb,        ""));
      AgentInfoVisitor(agentInfoWindow, _trackid).visit(item._pAgent[0]);
      AgentInfoVisitor(agentInfoWindow, _trackid).visit(item._pAgent[1]);

      for(size_t i=0; i<3; i++) //Hacky way to keep expanded information
        if(agentInfoWindow->topLevelItem(i))
          agentInfoWindow->topLevelItem(i)->setExpanded(topLevelExpanded[i]);
      agentInfoWindow->blockSignals(false);
      _pItfc->getSimulation().unlock();
    }
  else
    {
      agentInfoWindow->hide();
      _trackid = -1;
    }
  if (_selRow != -1 && _selCol != -1)
    {
      const ScalarAgentGrid* pAgentGrid = static_cast<ScalarAgentGrid*>(_pItfc->getScalarAgentGrid());
      const ScalarAgentItem& item
        = pAgentGrid->getGrid()[_selRow * _ui.glWidget->dim.x + _selCol];

      QString str = QString("pos = (%1,%2)").arg(_selCol).arg(_selRow);
#if 0
      if (GET_BIT(item._bitMask, ScalarAgentGrid::_bitCas))
        {
          str += ", caseous";
        }
      else
        {
          str += QString(", killings = %1").arg(item._nKillings);
        }

      if (GET_BIT(item._bitMask, ScalarAgentGrid::_bitSrc))
        {
          str += ", source";

          bool srcMac = GET_BIT(item._bitMask, ScalarAgentGrid::_bitSrcMac);
          bool srcTgam= GET_BIT(item._bitMask, ScalarAgentGrid::_bitSrcTgam);
          bool srcTcyt = GET_BIT(item._bitMask, ScalarAgentGrid::_bitSrcTcyt);
          bool srcTreg = GET_BIT(item._bitMask, ScalarAgentGrid::_bitSrcTreg);

          str += " (";

          bool first = false;
          if (srcMac)
            {
              str += "M";
              first = true;
            }

          if (srcTgam)
            {
              if (first) str += ",";
              else first = true;
              str += "Tg";
            }

          if (srcTcyt)
            {
              if (first) str += ",";
              else first = true;
              str += "Tc";
            }

          if (srcTreg)
            {
              if (first) str += ",";
              else first = true;
              str += "Tr";
            }
          str += ")";
        }

      str += QString(", Be = %1").arg(item._extMtb, 0, 'f', 2);
#endif

      glColor3f(1.0f, 0.0f, 0.0f);
      _ui.glWidget->renderText(10, 20, str);

      str = QString("(%1, %2, %3, %4, %5)").arg(item._attractant, 0, 'f', 2).
            arg(item._TNF, 0, 'f', 2).arg(item._CCL2, 0, 'f', 2).
            arg(item._CCL5, 0, 'f', 2).arg(item._CXCL9, 0, 'f', 2);
      _ui.glWidget->renderText(10, 40, str);

#if 0
      // For debugging granuloma boundary definition.
      GridCell cell = _pItfc->getSimulation().getGrGrid()(_selRow, _selCol);
      bool occupied = cell.isOccupied();
      str = QString("(%1, %2, %3").arg(item._TNF, 0, 'f', 2).
            arg(item._attractant, 0, 'f', 2).
            arg(item._extMtb, 0, 'f', 2);
      str += occupied ? ", O )" : ", NO)";

      _ui.glWidget->renderText(10, 60, str);
#endif

      str = getAgentStr(item._pAgent[0]);
      _ui.glWidget->renderText(10, _ui.glWidget->height() - 10, str);

      str = getAgentStr(item._pAgent[1]);
      _ui.glWidget->renderText(10, _ui.glWidget->height() - 30, str);
    }
}

QString GLWindow::getAgentStr(const Agent* pAgent)
{
  if(!pAgent) return "";
  QString res(QString::number(pAgent->getID())+" ");

  if (Mac::isMac(pAgent))
    {
      const Mac* pMac = dynamic_cast<const Mac*>(pAgent);
      switch (pMac->getState())
        {
        case Mac::MAC_RESTING:
          res += "Mr";
          break;
        case Mac::MAC_INFECTED:
          res += "Mi";
          break;
        case Mac::MAC_CINFECTED:
          res += "Mci";
          break;
        case Mac::MAC_ACTIVE:
        {
          res += "Ma";

          int days, hours, minutes;
          GrSimulation::convertSimTime(pMac->getActivationTime(), days, hours, minutes);
          res += QString(", act = %1d%2h%3m").arg(days, 3, 10, QChar('0')).
                 arg(hours, 2, 10, QChar('0')).arg(minutes, 2, 10, QChar('0'));

          break;
        }
        case Mac::MAC_DEAD:
          return res;
        }

      if (pMac->getNFkB())
        {
          res += ", NFkB";
        }

      if (pMac->getStat1())
        {
          res += ", stat1";
        }

      if (pMac->isDeactivated())
        {
          res += ", deact";
        }
      res += QString(", [sbTNF/sbIL10] = %1").arg(pMac->getsurfBoundTNFR1() / pMac->getsurfBoundIL10R(), 0, 'f', 2);
    }
  return res;
}

QImage GLWindow::grabFrameBuffer()
{
  return _ui.glWidget->renderPixmap().toImage();
}

void GLWindow::updateColorMapLabels(float min, float max)
{
  float f2 = (max - min) * .75 + min;
  float f3 = (max - min) * 0.5 + min;
  float f4 = (max - min) * .25 + min;

  _ui.labelLegend1->setText(QString::number(max, '.', _COLORMAP_LEGEND_PRECISION));
  _ui.labelLegend2->setText(QString::number(f2, '.', _COLORMAP_LEGEND_PRECISION));
  _ui.labelLegend3->setText(QString::number(f3, '.', _COLORMAP_LEGEND_PRECISION));
  _ui.labelLegend4->setText(QString::number(f4, '.', _COLORMAP_LEGEND_PRECISION));
  _ui.labelLegend5->setText(QString::number(min, '.', _COLORMAP_LEGEND_PRECISION));
}

void GLWindow::selectCell(int row, int col)
{
  if (_selRow == row && _selCol == col)
    {
      // unselect
      _selRow = _selCol = -1;
      agentInfoWindow->hide();
    }
  else
    {
      _selRow = row;
      _selCol = col;
      if(agentInfoWindow->isHidden()) agentInfoWindow->show();
    }
  updateSelectedCellStats();
  updateWindow();
}

void GLWindow::moveSelectionLeft()
{
  if (_selRow != -1 && _selCol != -1)
    {
      emit updateSelection(_selRow, (_selCol - 1 + _ui.glWidget->dim.x) % _ui.glWidget->dim.x);
    }
}

void GLWindow::moveSelectionRight()
{
  if (_selRow != -1 && _selCol != -1)
    {
      emit updateSelection(_selRow, (_selCol + 1 + _ui.glWidget->dim.x) % _ui.glWidget->dim.x);
    }
}

void GLWindow::moveSelectionUp()
{
  if (_selRow != -1 && _selCol != -1)
    {
      emit updateSelection((_selCol + 1 + _ui.glWidget->dim.y) % _ui.glWidget->dim.y, _selCol);
    }
}

void GLWindow::moveSelectionDown()
{
  if (_selRow != -1 && _selCol != -1)
    {
      emit updateSelection((_selCol - 1 + _ui.glWidget->dim.y) % _ui.glWidget->dim.y, _selCol);
    }
}

void GLWindow::visualize()
{
  _pItfc->visualize();
}

void GLWindow::printText()
{
  const Simulation* pSim = &(_pItfc->getSimulation());
  updateSelectedCellStats();

  glColor3f(1.0f, 0.0f, 0.0f);

  int xOffset = 200;

  if (_printTime)
    {
      _ui.glWidget->renderText(_ui.glWidget->width() - xOffset, 20,
                               Simulation::getTimeStr(pSim->getTime(), _pItfc->getTime()));
    }

  if (_printOutcome)
    {
      QString text;
      for (int i = 0; i < NOUTCOMES; i++)
        {
          if (i != 0)
            text += ", ";

          switch (_status[i])
            {
            case GR_NONE:
              text += "None";
              break;
            case GR_UNKNOWN:
              text += "Unknown";
              break;
            case GR_CLEARANCE:
              text += "Clearance";
              break;
            case GR_CONTAINMENT:
              text += "Containment";
              break;
            case GR_CONTAINMENT_INCONSISTENT:
              text += "Containment?";
              break;
            case GR_DISSEMINATION:
              text += "Dissemination";
              break;
            case GR_DISSEMINATION_INCONSISTENT:
              text += "Dissemination?";
              break;

            }
        }

      _ui.glWidget->renderText(_ui.glWidget->width() - xOffset, 40, text);
    }
}

void GLWindow::resizeGLWidget(int width, int height)
{
  int dWidth = this->width() - _ui.glWidget->width();
  int dHeight = this->height() - _ui.glWidget->height();

  resize(width + dWidth, height + dHeight);
}

void GLWindow::closeEvent(QCloseEvent* pEvent)
{
  pEvent->accept();
  qApp->quit();
}
