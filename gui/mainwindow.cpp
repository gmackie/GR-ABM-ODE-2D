#include "mainwindow.h"
#include "glwidget.h"
#include <assert.h>
#include <iostream>
#include <float.h>
#include <fstream>
#include <QFileDialog>
#include <QMessageBox>
#include <QShortcut>
#include <QSettings>
#include <vector>
#include <list>
#include <limits>

//TODO: Clean this up...

#include "simulation.h"
#include "colormaps/rainbow.h"
#include "colormaps/blackwhite.h"
#include "colormaps/fire.h"
#include "colormaps/coolwarm.h"
#include "colormaps/greenred.h"
#include "glyphs/glyph.h"
#include "glyphs/glyphhedgehog.h"
#include "glyphs/glyphcone.h"
#include "glyphs/glypharrow.h"
#include "glyphs/glyphtexture.h"
#include "scalardatasets/scalardataset.h"
#include "scalardatasets/scalarindexeddataset.h"
#include "scalardatasets/scalartnfattrextmtb.h"
#include "scalardatasets/scalarcelldensitydataset.h"
#include "scalardatasets/scalarintmtbdataset.h"
#include "scalardatasets/scalargrowthratedataset.h"
#include "scalardatasets/scalartotmtbdataset.h"
#include "scalardatasets/scalardivergencedataset.h"
#include "vectordatasets/vectorgradientdataset.h"
#include "vectordatasets/vector.h"
#include "gui/glwidget.h"
#include "gui/paramwindow.h"
#include "snapshot.h"
#include "colormaps/colormap.h"
#include "scalardatasets/scalarnormalizer.h"
#include "scalardatasets/scalargrid.h"
#include "gui/glwindow.h"
#include "gui/paramwindow.h"
#include "gui/statwidget.h"
#include "gui/agentswidget.h"
#include "gui/graphcontroller.h"

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

// static members
MainWindow* MainWindow::_pMainWindow = NULL;

const QString MainWindow::_INTERPOLATION_BILINEAR = "Bi-linear";
const QString MainWindow::_INTERPOLATION_NEAREST_NEIGHBOR = "Nearest neighbour";

const QString MainWindow::_GLYPH_HEDGEHOG = "Hedgehog";
const QString MainWindow::_GLYPH_CONE = "Cone";
const QString MainWindow::_GLYPH_ARROW = "Arrow";
const QString MainWindow::_GLYPH_TEXTURE = "Texture";

MainWindow::MainWindow(MainInterface* pItfc, GLWindow* pGLWindow, QWidget* pParamWindow,
                       StatWidget* pStatWidget, AgentsWidget* pAgentsWidget, const QDir& dir, QWidget* parent)
  : QMainWindow(parent)
  , _ui()
  , _pItfc(pItfc)
  , _pAgentHistogram(new AgentHistogram(pItfc->getSimulation()))
  , _pGraphController(new GraphController(pItfc->getSimulation().getStats(), 144*200))
  , _pGLWindow(pGLWindow)
  , _pParamWindow(pParamWindow)
  , _pStatWidget(pStatWidget)
  , _pAgentsWidget(pAgentsWidget)
  , _pCurrentColorMap(_pItfc->getColorMapSmoke())
  , _pCurrentNormalizer(_pItfc->getScalarSmokeNormalizer())
  , _pCurrentScalarGrid(_pItfc->getScalarSmokeGrid())
  , _timerId(0)
  , _simStatus(SIM_STOPPED)
  , _stopwatch(QTime::currentTime())
  , _pSnapshot(NULL)
  , _scriptingMode(false)
  , _dir(dir)
{
  assert(_pItfc);
  assert(_pGLWindow);
  assert(_pParamWindow);
  assert(_pStatWidget);
  assert(_pAgentsWidget);

  _pMainWindow = this;
  _ui.setupUi(this);
  _ui.tabStatistics->layout()->addWidget(_pStatWidget);
  //_ui.tabStatistics->layout()->addItem(new QSpacerItem(0, 0, QSizePolicy::Minimum, QSizePolicy::Expanding));
  _ui.tabAgents->layout()->addWidget(_pAgentsWidget);
  _ui.tabAgents->layout()->addItem(new QSpacerItem(0, 0, QSizePolicy::Minimum, QSizePolicy::Expanding));

  initSimulationTab();
  initVisualizationTab();
  initOutcomeTab();
  initOutputTab();
  initColorMapTab();
  initSmokeTab();
  initGlyphsTab();
  initIsolinesTab();
  initHeightTab();
  initShortcuts();
  loadSettings();

  _ui.labelTransparency->setEnabled(_ui.checkBoxBlend->isChecked());
  _ui.labelHeightGridTransparency->setEnabled(_ui.checkBoxBlend->isChecked());
  _ui.horizontalSliderAlpha->setEnabled(_ui.checkBoxBlend->isChecked());
  _ui.horizontalSliderHeightAlpha->setEnabled(_ui.checkBoxBlend->isChecked());
  _ui.tabSmoke->setEnabled(_ui.checkBoxDrawSmoke->isChecked());
  _ui.tabGlyphs->setEnabled(_ui.checkBoxDrawGlyphs->isChecked());
  _ui.tabHeightPlot->setEnabled(_ui.checkBoxDrawHeightPlot->isChecked());
  _ui.groupBoxIsoline->setEnabled(_ui.checkBoxDrawIsolines->isChecked());
  _ui.tabAgents->setEnabled(_ui.checkBoxDrawAgents->isChecked());

  /* set connections */
  connect(_ui.pushButtonAnimation, SIGNAL(clicked(bool)), this, SLOT(toggleAnimation(void)));
  connect(this, SIGNAL(updateGL(void)), _pGLWindow, SLOT(updateWindow(void)));
  connect(_pGLWindow, SIGNAL(updateSelection(int, int)), _pAgentsWidget, SLOT(setAgentSelection(int, int)));
  connect(_pGLWindow, SIGNAL(updateSelection(int, int)), _pGLWindow, SLOT(selectCell(int, int)));
  connect(this, SIGNAL(updateSelectedCellStats(void)), _pGLWindow, SLOT(updateSelectedCellStats(void)));
  connect(&_pItfc->getSimulation(), SIGNAL(stopConditionMet(void)), this, SLOT(stop(void)));
  connect(&_pItfc->getSimulation(), SIGNAL(updated(void)), this, SLOT(update(void)));
  connect(_pAgentsWidget, SIGNAL(updateGL(void)), _pGLWindow, SLOT(updateWindow(void)));
  connect(_ui.newGraphWindowButton, SIGNAL(clicked(void)), _pGraphController, SLOT(showNewTimeGraph(void)));
  connect(_pAgentHistogram, SIGNAL(closing(bool)), _ui.agentHistButton, SLOT(setChecked(bool)));
  connect(_pAgentsWidget, SIGNAL(agentFilterChanged(int, int, int, bool)), _pAgentHistogram, SLOT(setAgentFilter(int,int,int,bool)));

  int deltaHue = _ui.horizontalSliderHue->value();
  _pCurrentColorMap->setHueDelta(deltaHue / _SLIDER_COLORMAP_HUE);

  int deltaSat = _ui.horizontalSliderSat->value();
  _pCurrentColorMap->setSatDelta(deltaSat / _SLIDER_COLORMAP_SAT);

  int deltaVal = _ui.horizontalSliderVal->value();
  _pCurrentColorMap->setValDelta(deltaVal / _SLIDER_COLORMAP_VAL);

  int alpha = _ui.horizontalSliderAlpha->value();
  _pCurrentColorMap->setAlpha(alpha / _SLIDER_COLORMAP_ALPHA);

  int nColors = _ui.spinBoxNumberOfColors->value();
  _pCurrentColorMap->setNrBands(nColors);

  _pCurrentColorMap->setInvert(_ui.checkBoxReverseColorMap->isChecked());
  emit updateColorMap(_pCurrentColorMap);
}

void MainWindow::initOutcomeTab()
{

  for(int i=0;i<NOUTCOMES;i++)
    _ui.comboBoxOutcomeNumber->addItem(QString("Method %1").arg(i+1));

  for(int i=0;i<NOUTCOMEMETH;i++) {
    std::ostringstream ss;
    ss<<OutcomeMethod(i);
    _ui.comboBoxOutcomeMethod->addItem(QString::fromStdString(ss.str()));
  }

  connect(_ui.comboBoxOutcomeMethod, SIGNAL(currentIndexChanged(const QString&)), this, SLOT(updateOutcomeSettings(void)));
  connect(_ui.comboBoxOutcomeNumber, SIGNAL(currentIndexChanged(const QString&)), this, SLOT(updateOutcomeParameters(void)));
  connect(_ui.spinBoxOutcomeSamplePeriod, SIGNAL(valueChanged(int)), this, SLOT(updateOutcomeSettings(void)));
  connect(_ui.spinBoxOutcomeTestPeriod, SIGNAL(valueChanged(int)), this, SLOT(updateOutcomeSettings(void)));
  connect(_ui.doubleSpinBoxOutcomeAlpha, SIGNAL(valueChanged(double)), this, SLOT(updateOutcomeSettings(void)));
}

void MainWindow::initOutputTab()
{
  _ui.spinBoxSnapshotCsvInterval->setValue(_SNAPSHOT_CSV_PERIOD);
  _ui.spinBoxSnapshotPicInterval->setValue(_SNAPSHOT_PIC_PERIOD);

  connect(_ui.pushButtonDumpGrids, SIGNAL(clicked(bool)), this, SLOT(dumpGrids(void)));
  connect(_ui.pushButtonLoadState, SIGNAL(clicked(bool)), this, SLOT(loadState(void)));
  connect(_ui.pushButtonSaveState, SIGNAL(clicked(bool)), this, SLOT(saveState(void)));
  connect(_ui.pushButtonTakeSnapshot, SIGNAL(clicked(bool)), this, SLOT(takePictureSnapshot(void)));
}

void MainWindow::initHeightTab()
{
  /* configure _ui.comboBoxHeightDataset */
  initComboAllScalars(_ui.comboBoxHeightDataset);
  _ui.comboBoxHeightDataset->setCurrentIndex(GrGrid::IDX_TNF);

  /* configure _ui.comboBoxHeightColorDataset */
  initComboAllScalars(_ui.comboBoxHeightColorDataset);
  _ui.comboBoxHeightColorDataset->setCurrentIndex(GrGrid::IDX_TNF);

  connect(_ui.comboBoxHeightDataset, SIGNAL(currentIndexChanged(int)), this, SLOT(updateHeightPlotSettings(void)));
  connect(_ui.comboBoxHeightColorDataset, SIGNAL(currentIndexChanged(int)), this, SLOT(updateHeightPlotSettings(void)));
  connect(_ui.comboBoxHeightMappingMethod, SIGNAL(currentIndexChanged(int)), this, SLOT(updateHeightPlotSettings(void)));
  connect(_ui.pushButtonHeightRescale, SIGNAL(clicked(bool)), this, SLOT(rescaleHeightRequest(void)));
  connect(_ui.doubleSpinBoxHeightMin, SIGNAL(valueChanged(double)), this, SLOT(updateHeightMinMax(void)));
  connect(_ui.doubleSpinBoxHeightMax, SIGNAL(valueChanged(double)), this, SLOT(updateHeightMinMax(void)));
  connect(_ui.doubleSpinBoxHeightMaxHeight, SIGNAL(valueChanged(double)), this, SLOT(updateHeightPlotSettings(void)));
  connect(_ui.checkBoxHeightDrawGrid, SIGNAL(toggled(bool)), this, SLOT(updateHeightPlotSettings(void)));
  connect(_ui.doubleSpinBoxHeightMax_2, SIGNAL(valueChanged(double)), this, SLOT(updateHeightPlotSettings(void)));
  connect(_ui.horizontalSliderHeightAlpha, SIGNAL(valueChanged(int)), this, SLOT(updateHeightPlotSettings(void)));
}

void MainWindow::initShortcuts()
{
  QShortcut* pAnimate = new QShortcut(QKeySequence("space"), this, NULL, NULL, Qt::ApplicationShortcut);
  QShortcut* pAgents = new QShortcut(QKeySequence("a"), this, NULL, NULL, Qt::ApplicationShortcut);
  QShortcut* pSmoke = new QShortcut(QKeySequence("s"), this, NULL, NULL, Qt::ApplicationShortcut);
  QShortcut* pGlyphs = new QShortcut(QKeySequence("g"), this, NULL, NULL, Qt::ApplicationShortcut);
  QShortcut* pIsolines = new QShortcut(QKeySequence("i"), this, NULL, NULL, Qt::ApplicationShortcut);
  QShortcut* pHeightPlot = new QShortcut(QKeySequence("h"), this, NULL, NULL, Qt::ApplicationShortcut);
  QShortcut* pGranulomaBorder = new QShortcut(QKeySequence("b"), this, NULL, NULL, Qt::ApplicationShortcut);
  QShortcut* pSmokeColorMap = new QShortcut(QKeySequence("Ctrl+Shift+s"), this, NULL, NULL, Qt::ApplicationShortcut);
  QShortcut* pGlyphsColorMap = new QShortcut(QKeySequence("Ctrl+Shift+g"), this, NULL, NULL, Qt::ApplicationShortcut);
  QShortcut* pIsolinesColorMap = new QShortcut(QKeySequence("Ctrl+Shift+i"), this, NULL, NULL, Qt::ApplicationShortcut);
  QShortcut* pHeightPlotColorMap = new QShortcut(QKeySequence("Ctrl+Shift+h"), this, NULL, NULL, Qt::ApplicationShortcut);
  QShortcut* pColorMap0 = new QShortcut(QKeySequence("1"), this, NULL, NULL, Qt::ApplicationShortcut);
  QShortcut* pColorMap1 = new QShortcut(QKeySequence("2"), this, NULL, NULL, Qt::ApplicationShortcut);
  QShortcut* pColorMap2 = new QShortcut(QKeySequence("3"), this, NULL, NULL, Qt::ApplicationShortcut);
  QShortcut* pColorMap3 = new QShortcut(QKeySequence("4"), this, NULL, NULL, Qt::ApplicationShortcut);
  QShortcut* pColorMap4 = new QShortcut(QKeySequence("5"), this, NULL, NULL, Qt::ApplicationShortcut);
  QShortcut* pDumpGrids = new QShortcut(QKeySequence("d"), this, NULL, NULL, Qt::ApplicationShortcut);
  QShortcut* pReverseColorMap = new QShortcut(QKeySequence("r"), this, NULL, NULL, Qt::ApplicationShortcut);
  QShortcut* pQuit = new QShortcut(QKeySequence("q"), this, NULL, NULL, Qt::ApplicationShortcut);
  QShortcut* p2dView = new QShortcut(QKeySequence("F1"), this, NULL, NULL, Qt::ApplicationShortcut);
  QShortcut* p3dViewHeight = new QShortcut(QKeySequence("F3"), this, NULL, NULL, Qt::ApplicationShortcut);
  QShortcut* pBlend = new QShortcut(QKeySequence("Ctrl+b"), this, NULL, NULL, Qt::ApplicationShortcut);
  QShortcut* pArrowLeft = new QShortcut(QKeySequence("left"), this, NULL, NULL, Qt::ApplicationShortcut);
  QShortcut* pArrowRight = new QShortcut(QKeySequence("right"), this, NULL, NULL, Qt::ApplicationShortcut);
  QShortcut* pArrowUp = new QShortcut(QKeySequence("up"), this, NULL, NULL, Qt::ApplicationShortcut);
  QShortcut* pArrowDown = new QShortcut(QKeySequence("down"), this, NULL, NULL, Qt::ApplicationShortcut);
  QShortcut* pFullscreen = new QShortcut(QKeySequence("F11"), this, NULL, NULL, Qt::ApplicationShortcut);
  QShortcut* pLoadState = new QShortcut(QKeySequence("Ctrl+l"), this, NULL, NULL, Qt::ApplicationShortcut);
  QShortcut* pSaveState = new QShortcut(QKeySequence("Ctrl+s"), this, NULL, NULL, Qt::ApplicationShortcut);
  QShortcut* pTakePicture = new QShortcut(QKeySequence("F10"), this, NULL, NULL, Qt::ApplicationShortcut);

  connect(pAnimate, SIGNAL(activated(void)), _ui.pushButtonAnimation, SLOT(click(void)));
  connect(pAgents, SIGNAL(activated(void)), _ui.checkBoxDrawAgents, SLOT(click(void)));
  connect(pSmoke, SIGNAL(activated(void)), _ui.checkBoxDrawSmoke, SLOT(click(void)));
  connect(pGlyphs, SIGNAL(activated(void)), _ui.checkBoxDrawGlyphs, SLOT(click(void)));
  connect(pIsolines, SIGNAL(activated(void)), _ui.checkBoxDrawIsolines, SLOT(click(void)));
  connect(pHeightPlot, SIGNAL(activated(void)), _ui.checkBoxDrawHeightPlot, SLOT(click(void)));
  connect(pGranulomaBorder, SIGNAL(activated(void)), _ui.checkBoxDrawGranulomaBorder, SLOT(click(void)));
  connect(pSmokeColorMap, SIGNAL(activated(void)), this, SLOT(selectSmokeColorMap(void)));
  connect(pGlyphsColorMap, SIGNAL(activated(void)), this, SLOT(selectGlyphsColorMap(void)));
  connect(pIsolinesColorMap, SIGNAL(activated(void)), this, SLOT(selectIsolinesColorMap(void)));
  connect(pHeightPlotColorMap, SIGNAL(activated(void)), this, SLOT(selectHeightPlotColorMap(void)));
  connect(pColorMap0, SIGNAL(activated(void)), this, SLOT(selectColorMap0(void)));
  connect(pColorMap1, SIGNAL(activated(void)), this, SLOT(selectColorMap1(void)));
  connect(pColorMap2, SIGNAL(activated(void)), this, SLOT(selectColorMap2(void)));
  connect(pColorMap3, SIGNAL(activated(void)), this, SLOT(selectColorMap3(void)));
  connect(pColorMap4, SIGNAL(activated(void)), this, SLOT(selectColorMap4(void)));
  connect(pDumpGrids, SIGNAL(activated(void)), _ui.pushButtonDumpGrids, SLOT(click(void)));
  connect(pReverseColorMap, SIGNAL(activated(void)), _ui.checkBoxReverseColorMap, SLOT(click(void)));
  connect(pQuit, SIGNAL(activated(void)), qApp, SLOT(closeAllWindows(void)));
  connect(p2dView, SIGNAL(activated(void)), _pGLWindow, SIGNAL(set2DView(void)));
  connect(p3dViewHeight, SIGNAL(activated(void)), _pGLWindow, SIGNAL(set3DView(void)));
  connect(pBlend, SIGNAL(activated(void)), _ui.checkBoxBlend, SLOT(click(void)));
  connect(pArrowLeft, SIGNAL(activated(void)), _pGLWindow, SLOT(moveSelectionLeft(void)));
  connect(pArrowRight, SIGNAL(activated(void)), _pGLWindow, SLOT(moveSelectionRight(void)));
  connect(pArrowUp, SIGNAL(activated(void)), _pGLWindow, SLOT(moveSelectionUp(void)));
  connect(pArrowDown, SIGNAL(activated(void)), _pGLWindow, SLOT(moveSelectionDown(void)));
  connect(pFullscreen, SIGNAL(activated(void)), _pGLWindow, SLOT(toggleFullScreen(void)));
  connect(pLoadState, SIGNAL(activated(void)), this, SLOT(loadState(void)));
  connect(pSaveState, SIGNAL(activated(void)), this, SLOT(saveState(void)));
  connect(pTakePicture, SIGNAL(activated(void)), this, SLOT(takePictureSnapshot(void)));
}

void MainWindow::initVisualizationTab()
{
  _ui.checkBoxBlend->setChecked(_pItfc->getBlend());
  _ui.checkBoxDrawAgents->setChecked(_pItfc->getDrawAgents());
  _ui.checkBoxDrawGlyphs->setChecked(_pItfc->getDrawGlyphs());
  _ui.checkBoxDrawGranulomaBorder->setChecked(_pItfc->getDrawGranuloma());
  _ui.checkBoxDrawHeightPlot->setChecked(_pItfc->getDrawHeightPlot());
  _ui.checkBoxDrawIsolines->setChecked(_pItfc->getDrawIsolines());
  _ui.checkBoxDrawSmoke->setChecked(_pItfc->getDrawSmoke());
  _ui.checkBoxDrawTime->setChecked(_pGLWindow->getPrintTime());
  _ui.checkBoxDrawOutcomes->setChecked(_pGLWindow->getPrintOutcome());

  initComboAllScalars(_ui.comboBoxGranulomaDataset);

  // Set the initial granuloma dataset.
  int defaultidx = GrGrid::IDX_NGRIDS+3;  //TNF+ATTR+MTB
  _ui.comboBoxGranulomaDataset->setCurrentIndex(defaultidx);

  // Overrides the granuloma dataset defined in the MainInterface constructor.
  ScalarDataset* pScalarGranulomaDataset = getNewScalarDataset(defaultidx);
  _pItfc->setScalarGranulomaDataset(pScalarGranulomaDataset);

  connect(_ui.checkBoxDrawTime, SIGNAL(toggled(bool)), _pGLWindow, SLOT(setPrintTime(bool)));
  connect(_ui.checkBoxDrawOutcomes, SIGNAL(toggled(bool)), _pGLWindow, SLOT(setPrintOutcome(bool)));
  connect(_ui.checkBoxDrawGranulomaBorder, SIGNAL(toggled(bool)), this, SLOT(updateGranulomaSettings(void)));
  connect(_ui.comboBoxGranulomaDataset, SIGNAL(currentIndexChanged(int)), this, SLOT(updateGranulomaSettings(void)));
  connect(_ui.horizontalSliderGranulomaBorderThreshold, SIGNAL(valueChanged(int)), this, SLOT(updateGranulomaSettings(void)));

}

/**
 * @name Simulation Tab management
 * @{
 */

/**
 * @brief Initialize the simulation tab with default values
 */

void MainWindow::initSimulationTab()
{
  /* configure _ui.comboBoxDiffusion */
  // The first item added will cause a signal.
  // Block it til we assign the right diffusion method
  _ui.comboBoxDiffusion->blockSignals(true);
  for(int i=0;i<NDIFFUSEMETH;i++) {
    std::ostringstream ss;
    ss<<DiffusionMethod(i);
    _ui.comboBoxDiffusion->addItem(QString::fromStdString(ss.str()));
  }
  _ui.comboBoxDiffusion->blockSignals(false);

  // Set the default diffusion method to ADE swap.
  // This must be after connecting the currentIndexChanged signal
  // to the updateDiffusionMethod slot. Otherwise, setting
  // the combo box doesn't update the simulation object diffusion
  // method, so it remains whatever the default method is from
  // the simulation object constructor.
  _ui.comboBoxDiffusion->setCurrentIndex(DIFF_ADE_SWAP);

  // So seeds don't get truncated.
  _ui.spinBoxSeed->blockSignals(true);  //Block it in case of deserialization issues.
  _ui.spinBoxSeed->setMaximum(std::numeric_limits<int>::max());
  _ui.spinBoxSeed->setValue((int) g_Rand.getSeed());
  _ui.spinBoxSeed->blockSignals(false);
}


void MainWindow::on_checkBoxStopDays_toggled(bool checked)
{
  _pItfc->getSimulation().setTimeToSimulate(checked ? _ui.spinBoxStopTime->value() : -1);
}

void MainWindow::on_spinBoxStopTime_valueChanged(int arg1)
{
  _ui.spinBoxStopDays->blockSignals(true);  //Prevent double signaling
  _ui.spinBoxStopDays->setValue(arg1/TIME_STEPS_PER_DAY);
  _ui.spinBoxStopDays->blockSignals(false);
  _pItfc->getSimulation().setTimeToSimulate(arg1);
}

void MainWindow::on_spinBoxStopDays_valueChanged(int arg1)
{
  _ui.spinBoxStopTime->setValue(arg1 * TIME_STEPS_PER_DAY);
}

void MainWindow::on_checkBoxStopClearance_toggled(bool checked)
{
  _pItfc->getSimulation().setMtbClearance(checked);
}

void MainWindow::on_pushButtonShowParams_clicked()
{
  _pParamWindow->show();
}

void MainWindow::on_actionShowParameters_triggered()
{
  _pParamWindow->show();
}

void MainWindow::on_comboBoxDiffusion_currentIndexChanged(int index)
{
  _pItfc->getSimulation().setDiffusionMethod(DiffusionMethod(index));
}

void MainWindow::on_spinBoxDelay_valueChanged(int arg1)
{
  _pItfc->getSimulation().setDelay(arg1);
}

void MainWindow::on_spinBoxSeed_valueChanged(int arg1)
{
  // Need to make sure the simulation thread doesn't touch the RNG while changing the seed.
  _pItfc->getSimulation().lock();
  _pItfc->getSimulation().modelLock();
    g_Rand.setSeed(arg1);
  _pItfc->getSimulation().modelUnlock();
  _pItfc->getSimulation().unlock();
}

///@}

void MainWindow::initColorMapTab()
{
  /* configure _ui.comboBoxScalarColoring */
  for(int i=0;i<NCOLORMAPS;i++) {
    std::ostringstream ss;
    ss<<COLORMAP(i);
    _ui.comboBoxScalarColoring->addItem(QString::fromStdString(ss.str()));
  }

  /* configure _ui.comboBoxColorMapping */
  _ui.comboBoxColorMapSource->blockSignals(true);
  for(int i=0;i<NCOLORMAPSRCS;i++) {
      std::ostringstream ss;
      ss<<COLORMAP_SRC(i);
      _ui.comboBoxColorMapSource->addItem(QString::fromStdString(ss.str()));
  }
  _ui.comboBoxColorMapSource->blockSignals(false);
  _ui.doubleSpinBoxMin->setDecimals(std::numeric_limits<double>::digits10);
  _ui.doubleSpinBoxMax->setDecimals(std::numeric_limits<double>::digits10);
  connect(_ui.comboBoxScalarColoring, SIGNAL(currentIndexChanged(int)), this, SLOT(setColorMap(int)));
  connect(_ui.spinBoxNumberOfColors, SIGNAL(valueChanged(int)), this, SLOT(updateColorMap(void)));
  connect(_ui.horizontalSliderHue, SIGNAL(valueChanged(int)), this, SLOT(updateColorMap(void)));
  connect(_ui.horizontalSliderSat, SIGNAL(valueChanged(int)), this, SLOT(updateColorMap(void)));
  connect(_ui.horizontalSliderVal, SIGNAL(valueChanged(int)), this, SLOT(updateColorMap(void)));
  connect(_ui.horizontalSliderAlpha, SIGNAL(valueChanged(int)), this, SLOT(updateColorMap(void)));
  connect(_ui.checkBoxReverseColorMap, SIGNAL(toggled(bool)), this, SLOT(updateColorMap(void)));
  connect(_ui.comboBoxMappingMethod, SIGNAL(currentIndexChanged(int)), this, SLOT(setMappingMethod(int)));
  connect(_ui.doubleSpinBoxMin, SIGNAL(valueChanged(double)), this, SLOT(updateColorMapLabels(void)));
  connect(_ui.doubleSpinBoxMax, SIGNAL(valueChanged(double)), this, SLOT(updateColorMapLabels(void)));
  connect(_ui.pushButtonRescale, SIGNAL(clicked(bool)), this, SLOT(rescaleRequest(void)));
  connect(this, SIGNAL(updateColorMap(ColorMap*)), _pGLWindow, SLOT(setColorMap(ColorMap*)));
  connect(this, SIGNAL(updateColorMapLabels(float, float)), _pGLWindow, SLOT(updateColorMapLabels(float, float)));
}

void MainWindow::initSmokeTab()
{
  /* configure _ui.comboBoxSmokeDataset */
  initComboAllScalars(_ui.comboBoxSmokeDataset);
  _ui.comboBoxSmokeDataset->setCurrentIndex(GrGrid::IDX_TNF);
}

void MainWindow::initGlyphsTab()
{
  /* configure _ui.comboBoxGlyphScalarDataset */
  initComboAllScalars(_ui.comboBoxGlyphScalarDataset);
  _ui.comboBoxGlyphScalarDataset->setCurrentIndex(GrGrid::IDX_TNF);


  /* configure _ui.comboBoxGlyphVectorDataset */
  initComboAllScalars(_ui.comboBoxGlyphVectorDataset);
  for(int i=0;i<_ui.comboBoxGlyphVectorDataset->count();i++)  //Add ' grad' to the end of the label (consistency)
    _ui.comboBoxGlyphVectorDataset->setItemText(i, _ui.comboBoxGlyphVectorDataset->itemText(i) + " grad");
  _ui.comboBoxGlyphVectorDataset->setCurrentIndex(GrGrid::IDX_TNF);


  /* configure _ui.comboBoxGlyphType */
  _ui.comboBoxGlyphType->addItem(_GLYPH_HEDGEHOG);
  _ui.comboBoxGlyphType->addItem(_GLYPH_CONE);
  _ui.comboBoxGlyphType->addItem(_GLYPH_ARROW);
  _ui.comboBoxGlyphType->addItem(_GLYPH_TEXTURE);

  _ui.doubleSpinBoxGlyphScaling->setMaximum(std::numeric_limits<double>::max());
  connect(_ui.comboBoxGlyphType, SIGNAL(currentIndexChanged(const QString&)), this, SLOT(setGlyph(const QString&)));
}

void MainWindow::initIsolinesTab()
{
  connect(_ui.horizontalSliderIsolinesMin, SIGNAL(valueChanged(int)), this, SLOT(updateIsolinesSettings(void)));
  connect(_ui.horizontalSliderIsolinesMax, SIGNAL(valueChanged(int)), this, SLOT(updateIsolinesSettings(void)));
  connect(_ui.spinBoxNumberOfIsolines, SIGNAL(valueChanged(int)), this, SLOT(updateIsolinesSettings(void)));
  connect(_ui.doubleSpinBoxIsolinesLineWidth, SIGNAL(valueChanged(double)), this, SLOT(updateIsolinesSettings(void)));
  connect(_ui.doubleSpinBoxMin, SIGNAL(valueChanged(double)), this, SLOT(updateIsolinesSettings(void)));
  connect(_ui.doubleSpinBoxMax, SIGNAL(valueChanged(double)), this, SLOT(updateIsolinesSettings(void)));
}

/*
	Initialize a combo box with all the scalar data sets.
*/
void MainWindow::initComboAllScalars(QComboBox *&comboBox)
{
  for (int i = 0; i < GrGrid::IDX_NGRIDS; i++)
    {
      std::ostringstream ss;
      ss<<GrGrid::GRID_IDX(i);
      comboBox->addItem(QString::fromStdString(ss.str()));
    }
  comboBox->addItem("Cell Density");
  comboBox->addItem("Int. Mtb");
  comboBox->addItem("Tot. Mtb");
  comboBox->addItem("TNF + Attract + Ext. Mtb");
  comboBox->addItem("Internal Growth Rate");
}

void MainWindow::getScalarDataSetNames(std::string &names)
{
  std::ostringstream ss;
  for (int i = 0; i < GrGrid::IDX_NGRIDS; i++)
      ss<<GrGrid::GRID_IDX(i)<<std::endl;
  ss<<("Cell Density")<<std::endl;
  ss<<("Ext. Mtb")<<std::endl;
  ss<<("Int. Mtb")<<std::endl;
  ss<<("Tot. Mtb")<<std::endl;
  ss<<"TNF + Attract + Ext. Mtb"<<std::endl;
  ss<<"Internal Growth Rate"<<std::endl;
  names = ss.str();
}

MainWindow::~MainWindow()
{
  delete _pGraphController;
  delete _pAgentHistogram;
  delete _pSnapshot;
  _pMainWindow = NULL;

  Simulation& sim = _pItfc->getSimulation();

  // stop the simulation
  if (sim.isRunning())
    {
      // unlock the simulation thread if it is locked due to pausing
      if (_simStatus != SIM_RUNNING)
        sim.unlock();

      // signal simulation thread to stop
      sim.quit();

      // wait until the thread is stopped
      sim.wait();
    }
}

void MainWindow::update()
{
  timerEvent(NULL);
}

void MainWindow::timerEvent(QTimerEvent*)
{
  Simulation& sim = _pItfc->getSimulation();
  if(!sim.getUpdated()) return; // Nothing's changed, no need to update simulation data structures

  /* BEGIN CRITICAL SECTION */
  /* Since the mutex is recursive, a thread can lock the same mutex multiple times
   * and the mutex won't be unlocked until a corresponding number of unlock() calls have been made. */
  sim.lock();
  sim.setUpdated(false);
  _pItfc->doStep();

  int simTime = sim.getTime();
  _pItfc->setSimTime(simTime);
  const Stats stats = sim.getStats();

  //if(simTime % 72 == 0)   //Update every half day
  _pGraphController->update(stats, 1.0*simTime);
  if(!_pAgentHistogram->isHidden())
    _pAgentHistogram->updatePlot();

  _pItfc->incTime(_stopwatch.restart());

  updateWindowStatus();

  _pStatWidget->updateLabels(stats);
  emit updateSelectedCellStats();

  if (!_pCurrentNormalizer->getClamping())
    {
      float min = _pCurrentNormalizer->getMin();
      float max = _pCurrentNormalizer->getMax();

      _ui.doubleSpinBoxMin->setValue(min);
      _ui.doubleSpinBoxMax->setValue(max);
    }

  if (_pItfc->getDrawHeightPlot() && !_pItfc->getScalarHeightNormalizer()->getClamping())
    {
      float heightMin = _pItfc->getScalarHeightNormalizer()->getMin();
      float heightMax = _pItfc->getScalarHeightNormalizer()->getMax();

      _ui.doubleSpinBoxHeightMin->setValue(heightMin);
      _ui.doubleSpinBoxHeightMax->setValue(heightMax);
    }
  if (_pItfc->getDrawIsolines() && !_pItfc->getScalarIsolinesNormalizer()->getClamping())
    {
      int min = _ui.horizontalSliderIsolinesMin->value();
      int max = _ui.horizontalSliderIsolinesMax->value();

      float denormalizedMin = _pItfc->getScalarIsolinesNormalizer()->denormalize(min / 20.0f);
      float denormalizedMax = _pItfc->getScalarIsolinesNormalizer()->denormalize(max / 20.0f);
      _ui.labelIsolinesMin->setText(QString::number(denormalizedMin, '.', 4));
      _ui.labelIsolinesMax->setText(QString::number(denormalizedMax, '.', 4));
    }

  emit updateGL();

  // Take snapshots if requested and it is time.
  if (_pSnapshot)
    {
      // So the simulation doesn't advance while we are saving data.
      // Specifically needed when saving the simulation state, since that
      // uses the GrSimulation object, not the data saved from GrSimulation
      // in the MainInterface object.
      sim.modelLock();

      int picInterval = _ui.spinBoxSnapshotPicInterval->value();
      int csvInterval = _ui.spinBoxSnapshotCsvInterval->value();
      int stateInterval = _ui.spinBoxSnapshotStateInterval->value();
      simTime = sim.getTime();  //Get time again since this could have changed

      if (isTimeForAction(picInterval, simTime))
        {
          _pSnapshot->takePicture(simTime, _pGLWindow->grabFrameBuffer());
        }

      if (isTimeForAction(csvInterval, simTime))
        {
          _pSnapshot->takeSnapshot(simTime, stats);
        }

      if (isTimeForAction(stateInterval, simTime))
        {
          _pSnapshot->takeStateSnapshot(simTime, sim);
        }

      sim.modelUnlock();
    }
  sim.unlock();

}

// Determine whether or not an action should be taken on this time step,
// based on an action interval. Actions are generally not part of the
// model calculation, but doing something with the model state after
// it has been calculated for a time step, such as saving statistics
// or saving a graphics snapshot, etc.
//
// If the action interval is 0 the action is disabled.
// Otherwise if that interval of time steps has elapsed since the
// last time the action was taken then take the action again.
// We also do the action on the last time step of the simulation,
// whether that last time step occurred because the simulation reached
// its time limit or is finished for some other reason, such as
// clearance occurred.

bool MainWindow::isTimeForAction(int actionInterval, int simTime)
{
  bool result = (actionInterval > 0 && ((simTime % actionInterval == 0) || _simStatus == SIM_STOPPED));
  return result;
}

void MainWindow::updateWindowStatus()
{
  Simulation& sim = _pItfc->getSimulation();
  const Stats stats = sim.getStats();

  for (int i = 0; i < NOUTCOMES; i++)
    {
      _pGLWindow->setGrStatus(i, stats.getStatus(i));
    }
}

void MainWindow::toggleAnimation()
{
  if (_simStatus == SIM_RUNNING)
    switchStatus(SIM_PAUSED);
  else if (_simStatus == SIM_PAUSED)
    switchStatus(SIM_RUNNING);
  else
    switchStatus(SIM_RUNNING);
}

void MainWindow::updateColorMap()
{
  setColorMap(_ui.comboBoxScalarColoring->currentIndex());
}

void MainWindow::updateColorMapLabels()
{
  float min = _ui.doubleSpinBoxMin->value();
  float max = _ui.doubleSpinBoxMax->value();

  _pCurrentNormalizer->setMin(min);
  _pCurrentNormalizer->setMax(max);

  emit updateColorMapLabels(min, max);
  emit updateGL();
}

void MainWindow::rescaleRequest()
{
  float min = _pCurrentScalarGrid->getMin();
  float max = _pCurrentScalarGrid->getMax();

  _pCurrentNormalizer->setMin(min);
  _pCurrentNormalizer->setMax(max);

  _ui.doubleSpinBoxMin->setValue(min);
  _ui.doubleSpinBoxMax->setValue(max);

  emit updateGL();
}

void MainWindow::rescaleHeightRequest()
{
  float heightMin = _pItfc->getScalarHeightGrid()->getMin();
  float heightMax = _pItfc->getScalarHeightGrid()->getMax();

  _ui.doubleSpinBoxHeightMin->setValue(heightMin);
  _ui.doubleSpinBoxHeightMax->setValue(heightMax);

  emit updateGL();
}

void MainWindow::setMappingMethod(int value)
{
  if (value == 0)
    {
      _ui.doubleSpinBoxMax->setEnabled(true);
      _ui.doubleSpinBoxMax->setDecimals(3);
      _ui.doubleSpinBoxMin->setEnabled(true);
      _ui.doubleSpinBoxMin->setDecimals(3);
      _ui.pushButtonRescale->setEnabled(true);

      _pCurrentNormalizer->setClamping(true);

      // no need to update, as the previous mode was _MAP_METHOD_SCALE
    }
  else if (value == 1)
    {
      _ui.doubleSpinBoxMax->setEnabled(false);
      _ui.doubleSpinBoxMax->setDecimals(7);
      _ui.doubleSpinBoxMin->setEnabled(false);
      _ui.doubleSpinBoxMin->setDecimals(7);
      _ui.pushButtonRescale->setEnabled(false);

      _pCurrentNormalizer->setClamping(false);
      rescaleRequest();
    }
}

void MainWindow::setColorMap(const int value)
{
  switch(value) {
    case CMAPRAINBOW:
      _pCurrentColorMap = new ColorMapRainbow(); break;
    case CMAPGRAY:
      _pCurrentColorMap = new ColorMapBlackWhite(); break;
    case CMAPFIRE:
      _pCurrentColorMap = new ColorMapFire(); break;
    case CMAPCOOLWARM:
      _pCurrentColorMap = new ColorMapCoolWarm(); break;
    case CMAPGREENRED:
      _pCurrentColorMap = new ColorMapGreenRed(); break;
    default:
      assert(!"Invalid color map"); break;
    }
  assert(_pCurrentColorMap);

  const float deltaHue = 1.0f*_ui.horizontalSliderHue->value();
  _pCurrentColorMap->setHueDelta(deltaHue / _ui.horizontalSliderHue->maximum());

  const float deltaSat = 1.0f*_ui.horizontalSliderSat->value();
  _pCurrentColorMap->setSatDelta(deltaSat / _ui.horizontalSliderSat->maximum());

  const float deltaVal = 1.0f*_ui.horizontalSliderVal->value();
  _pCurrentColorMap->setValDelta(deltaVal / _ui.horizontalSliderVal->maximum());

  const float alpha = 1.0f*_ui.horizontalSliderAlpha->value();
  _pCurrentColorMap->setAlpha(alpha / _ui.horizontalSliderAlpha->maximum());

  const int nColors = _ui.spinBoxNumberOfColors->value();
  _pCurrentColorMap->setNrBands(nColors);

  _pCurrentColorMap->setInvert(_ui.checkBoxReverseColorMap->isChecked());

  const int colorMapSelection = _ui.comboBoxColorMapSource->currentIndex();
  switch(colorMapSelection) {
    case SMOKE_SRC:      _pItfc->setColorMapSmoke(_pCurrentColorMap); break;
    case GLYPHS_SRC:     _pItfc->setColorMapGlyphs(_pCurrentColorMap); break;
    case ISOLINES_SRC:   _pItfc->setColorMapIsolines(_pCurrentColorMap); break;
    case HEIGHTPLOT_SRC: _pItfc->setColorMapHeightPlot(_pCurrentColorMap); break;
    default: assert(!"Invalid colormap source"); break;
  }

  emit updateColorMap(_pCurrentColorMap);
  emit updateGL();
}

void MainWindow::setGlyph(const QString& value)
{
  Glyph* pGlyph = NULL;

  if (value == _GLYPH_HEDGEHOG)
    {
      pGlyph = new GlyphHedgehog();
    }
  else if (value == _GLYPH_CONE)
    {
      pGlyph = new GlyphCone();
    }
  else if (value == _GLYPH_ARROW)
    {
      pGlyph = new GlyphArrow();
    }
  else if (value == _GLYPH_TEXTURE)
    {
      QString fileName = QFileDialog::getOpenFileName(this, "Open Texture File", "", "");
      if (fileName == QString::null)
        {
          return _ui.comboBoxGlyphType->setCurrentIndex(0);
        }
      else
        {
          pGlyph = new GlyphTexture(fileName);
        }
    }
  else
    {
      assert(false);
    }

  _pItfc->setGlyph(pGlyph);

  emit updateGL();
}

void MainWindow::dumpGrids()
{
  QString fileName = QFileDialog::getSaveFileName(this, "Save Grids", "", "*.csv");
  if (fileName != QString::null)
    {
      const char* strFileName = fileName.toLatin1().data();
      std::cout << strFileName << std::endl;
      std::ofstream outFile(strFileName, std::ios_base::trunc);
      if (outFile)
        {
          _pItfc->serializeGrids(outFile);
          outFile.flush();
          outFile.close();
        }
      else
        {
          QMessageBox::critical(this, "Lung ABM",
                                "Failed to open/create file '" + fileName + "' for writing.\n"
                                "Please provide a non-existing file name.",
                                QMessageBox::Ok, QMessageBox::Ok);
        }
    }
}

void MainWindow::updateIsolinesSettings()
{
  float lineWidth = _ui.doubleSpinBoxIsolinesLineWidth->value();
  int nrIsolines = _ui.spinBoxNumberOfIsolines->value();

  int min = _ui.horizontalSliderIsolinesMin->value();
  int max = _ui.horizontalSliderIsolinesMax->value();

  float denormalizedMin = _pItfc->getScalarIsolinesNormalizer()->denormalize(min / _SLIDER_ISOLINES_MIN);
  float denormalizedMax = _pItfc->getScalarIsolinesNormalizer()->denormalize(max / _SLIDER_ISOLINES_MAX);
  _ui.labelIsolinesMin->setText(QString::number(denormalizedMin, '.', _SLIDER_ISOLINES_LABEL_PRECISION));
  _ui.labelIsolinesMax->setText(QString::number(denormalizedMax, '.', _SLIDER_ISOLINES_LABEL_PRECISION));

  _pItfc->setIsolinesWidth(lineWidth);

  if (nrIsolines == 1)
    {
      _ui.horizontalSliderIsolinesMax->setEnabled(false);
      _ui.horizontalSliderIsolinesMax->setValue(min);

      std::vector<float> vector(1, min / _SLIDER_ISOLINES_MIN);
      _pItfc->setIsolinesValueVector(vector);
    }
  else
    {
      _ui.horizontalSliderIsolinesMax->setEnabled(true);

      float delta = (float) (max/_SLIDER_ISOLINES_MAX - min/_SLIDER_ISOLINES_MIN) / (float) (nrIsolines - 1);

      std::vector<float> vector(nrIsolines, min / _SLIDER_ISOLINES_MIN);
      for (int i = 0; i < nrIsolines; i++)
        {
          vector[i] += (delta * i);
        }
      _pItfc->setIsolinesValueVector(vector);
    }

  emit updateGL();
}

void MainWindow::updateHeightPlotSettings()
{
  ScalarDataset* pScalarHeightDataset = getNewScalarDataset(_ui.comboBoxHeightDataset->currentIndex());
  ScalarDataset* pScalarHeightColorDataset = getNewScalarDataset(_ui.comboBoxHeightColorDataset->currentIndex());
  _pItfc->setScalarHeightDataset(pScalarHeightDataset);
  _pItfc->setScalarHeightColorDataset(pScalarHeightColorDataset);

  if (_ui.comboBoxHeightMappingMethod->currentIndex() == 0)
    {
      _ui.doubleSpinBoxHeightMax->setEnabled(true);
      _ui.doubleSpinBoxHeightMax->setDecimals(_SLIDER_HEIGHTPLOT_LABEL_PRECISION_CLAMPING);
      _ui.doubleSpinBoxHeightMin->setEnabled(true);
      _ui.doubleSpinBoxHeightMin->setDecimals(_SLIDER_HEIGHTPLOT_LABEL_PRECISION_CLAMPING);
      _ui.pushButtonHeightRescale->setEnabled(true);

      _pItfc->getScalarHeightNormalizer()->setClamping(true);
      _pItfc->getScalarHeightNormalizer()->setMin(_ui.doubleSpinBoxHeightMin->value());
      _pItfc->getScalarHeightNormalizer()->setMax(_ui.doubleSpinBoxHeightMax->value());
    }
  else if(_ui.comboBoxHeightMappingMethod->currentIndex() == 1)
    {
      _ui.doubleSpinBoxHeightMax->setEnabled(false);
      _ui.doubleSpinBoxHeightMax->setDecimals(_SLIDER_HEIGHTPLOT_LABEL_PRECISION_SCALING);
      _ui.doubleSpinBoxHeightMin->setEnabled(false);
      _ui.doubleSpinBoxHeightMin->setDecimals(_SLIDER_HEIGHTPLOT_LABEL_PRECISION_SCALING);
      _ui.pushButtonHeightRescale->setEnabled(false);

      _pItfc->getScalarHeightNormalizer()->setClamping(false);
      rescaleHeightRequest();
    }

  _pItfc->setHeightMax(_ui.doubleSpinBoxHeightMaxHeight->value());

  _pItfc->setDrawHeightGrid(_ui.checkBoxHeightDrawGrid->isChecked());
  _pItfc->setGridHeightMax(_ui.doubleSpinBoxHeightMax_2->value());
  _pItfc->setHeightGridAlpha(_ui.horizontalSliderHeightAlpha->value() / _SLIDER_HEIGHTPLOT_ALPHA);

  _pItfc->updateGrids();
  emit updateGL();
}

ScalarDataset* MainWindow::getNewScalarDataset(int idx)
{
  ScalarDataset* pScalarDataset = NULL;
  switch(idx) {
    case GrGrid::IDX_NGRIDS: pScalarDataset = new ScalarCellDensityDataset(); break;
    case GrGrid::IDX_NGRIDS+1: pScalarDataset = new ScalarIntMtbDataset();    break;
    case GrGrid::IDX_NGRIDS+2: pScalarDataset = new ScalarTotMtbDataset();    break;
    case GrGrid::IDX_NGRIDS+3: pScalarDataset = new ScalarTnfAttrExtMtb();    break;
    case GrGrid::IDX_NGRIDS+4: pScalarDataset = new ScalarGrowthRateDataset();break;
    default:
      if(idx < GrGrid::IDX_NGRIDS && idx >= 0)
        pScalarDataset = new ScalarIndexedDataset(GrGrid::GRID_IDX(idx));
      else
        assert(!"Invalid Scalar dataset index");
      break;
  }

  return pScalarDataset;
}

VectorDataset* MainWindow::getNewVectorDataset(int idx)
{
  return new VectorGradientDataset(getNewScalarDataset(idx), _pItfc->getVectorGlyphGrid());
}

void MainWindow::updateHeightMinMax()
{
  float min = _ui.doubleSpinBoxHeightMin->value();
  float max = _ui.doubleSpinBoxHeightMax->value();

  _pItfc->getScalarHeightNormalizer()->setMin(min);
  _pItfc->getScalarHeightNormalizer()->setMax(max);

  emit updateGL();
}

void MainWindow::stop()
{
  Simulation& sim = _pItfc->getSimulation();

  sim.lock();

  _pItfc->doStep();
  _pItfc->setSimTime(sim.getTime());
  _pStatWidget->updateLabels(sim.getStats());

  sim.unlock();

  switchStatus(SIM_STOPPED);

  _pItfc->incTime(_stopwatch.restart());

  updateWindowStatus();
  emit updateGL();

  if (_scriptingMode)
    qApp->quit();
}

void MainWindow::updateGranulomaSettings()
{
  _pItfc->setDrawGranuloma(_ui.checkBoxDrawGranulomaBorder->isChecked());

  ScalarDataset* pScalarGranulomaDataset = getNewScalarDataset(_ui.comboBoxGranulomaDataset->currentIndex());
  _pItfc->setScalarGranulomaDataset(pScalarGranulomaDataset);

  float val = _ui.horizontalSliderGranulomaBorderThreshold->value() / _SLIDER_GRANULOMA_BORDER_THRESHOLD;
  _pItfc->setGranulomaThreshold(val);

  emit updateGL();
}

void MainWindow::closeEvent(QCloseEvent*)
{
  if (_pGLWindow->isWindow())
    _pGLWindow->close();
}

void MainWindow::updateOutcomeParameters()
{
  Simulation& sim = _pItfc->getSimulation();

  int index = _ui.comboBoxOutcomeNumber->currentIndex();

  double alpha;
  int testPeriod, samplePeriod;

  OutcomeMethod method = sim.getOutcomeMethod(index);
  sim.getOutcomeParameters(index, samplePeriod, testPeriod, alpha);

  _ui.spinBoxOutcomeSamplePeriod->setValue(samplePeriod);
  _ui.spinBoxOutcomeTestPeriod->setValue(testPeriod);
  _ui.doubleSpinBoxOutcomeAlpha->setValue(alpha);

  bool checked = method != OUTCOME_NONE;
  _ui.labelOutcomeSamplePeriod->setEnabled(checked);
  _ui.labelOutcomeTestPeriod->setEnabled(checked);
  _ui.labelOutcomeAlpha->setEnabled(checked);
  _ui.spinBoxOutcomeSamplePeriod->setEnabled(checked);
  _ui.spinBoxOutcomeTestPeriod->setEnabled(checked);
  _ui.doubleSpinBoxOutcomeAlpha->setEnabled(checked);

  _ui.comboBoxOutcomeMethod->setCurrentIndex(method);
}

void MainWindow::updateOutcomeSettings()
{
  Simulation& sim = _pItfc->getSimulation();

  int index = _ui.comboBoxOutcomeNumber->currentIndex();

  OutcomeMethod method = OutcomeMethod(_ui.comboBoxOutcomeMethod->currentIndex());

  bool checked = method != OUTCOME_NONE;
  _ui.labelOutcomeSamplePeriod->setEnabled(checked);
  _ui.labelOutcomeTestPeriod->setEnabled(checked);
  _ui.labelOutcomeAlpha->setEnabled(checked);
  _ui.spinBoxOutcomeSamplePeriod->setEnabled(checked);
  _ui.spinBoxOutcomeTestPeriod->setEnabled(checked);
  _ui.doubleSpinBoxOutcomeAlpha->setEnabled(checked);

  int samplePeriod = _ui.spinBoxOutcomeSamplePeriod->value();
  int testPeriod = _ui.spinBoxOutcomeTestPeriod->value();
  double alpha = _ui.doubleSpinBoxOutcomeAlpha->value();

  sim.setOutcomeMethod(index, method, alpha, testPeriod, samplePeriod);
}

void MainWindow::loadState()
{
  QString fileName = QFileDialog::getOpenFileName(this, "Load state", _dir.absolutePath(), "Saved State files (*.state *.state.gz)");
  if (fileName != QString::null)
    {
      loadState(fileName.toLatin1().data());
    }
}


void MainWindow::loadState(std::string fileName)
{
  Simulation& sim = _pItfc->getSimulation();

  // pause the simulation
  if (_simStatus == SIM_RUNNING)
    switchStatus(SIM_PAUSED);

  std::ifstream in(fileName.c_str());
  if (!in.good())
    {
      std::string errstr = "Failed to open file '" + fileName + "'.";
      QString errmsg(errstr.c_str());
      QMessageBox::critical(this, "Lung ABM", errmsg, QMessageBox::Ok, QMessageBox::Ok);
    }
  else
    {
      boost::iostreams::filtering_istream in;

      QString qFileName = fileName.c_str();

      if(QFileInfo(qFileName).suffix() == "gz")
        {
          in.push(boost::iostreams::gzip_decompressor());
        }
      in.push(boost::iostreams::file_source(qFileName.toStdString()));

      sim.loadState(in);

      // update diffusion and seed in the user interface
      // The seed in the saved state (and in the Rand class) is of type unsigned int.
      // In _ui.spinBoxSeed, which is a QSpinBox obejct, it is int. So seeds from a saved
      // state might be truncated, most likely if they were originally based on a time value,
      // which is typically unsigned int. This is won't be a problem until about the year 2037.
      // The max 32 bit int is 2147483647. The number of seconds since 1970 (how unix time values
      // are defined) doesn't reach 2147483647 until about the year 2037.
      _ui.spinBoxSeed->blockSignals(true);
      _ui.spinBoxSeed->setValue((int) g_Rand.getSeed());
      _ui.spinBoxSeed->blockSignals(false);
      _ui.comboBoxDiffusion->blockSignals(true);
      _ui.comboBoxDiffusion->setCurrentIndex(sim.getDiffusionMethod());
      _ui.comboBoxDiffusion->blockSignals(false);

      // process/visualize new state
      timerEvent(NULL);
    }
}

void MainWindow::saveState()
{
  Simulation& sim = _pItfc->getSimulation();

  // pause the simulation, so that we don't need to lock
  if (_simStatus == SIM_RUNNING)
    switchStatus(SIM_PAUSED);

  QString fileName = QFileDialog::getSaveFileName(this, "Save state", _dir.absolutePath(), "Saved State files (*.state *.state.gz)");
  if (fileName != QString::null)
    {

      fileName += ".state.gz";
      namespace bio = boost::iostreams;
      bio::filtering_ostream out;
      out.push(bio::gzip_compressor());
      out.push(bio::file_sink(fileName.toStdString()));

      if (!out.good())
        {
          QMessageBox::critical(this, "Lung ABM",
                                "Failed to save file '" + fileName + "'.",
                                QMessageBox::Ok, QMessageBox::Ok);
        }
      else
        {
          sim.saveState(out);
          out.flush();
        }
    }
}

void MainWindow::takePictureSnapshot()
{
  QString fileName = QFileDialog::getSaveFileName(this, "Save picture", "", "*.png");
  if (fileName != QString::null)
    {
      QImage image = _pGLWindow->grabFrameBuffer();
      image.save(fileName, "png");
    }
}

void MainWindow::switchStatus(SimStatus newStatus)
{
  Simulation& sim = _pItfc->getSimulation();

  if(newStatus == _simStatus)
    return;

  switch(newStatus)
  {
    case SIM_RUNNING:
      if(sim.isRunning())
        sim.unlock();
      else
        sim.start();
      _pGLWindow->setWindowTitle("Visualization (running)");
      _ui.pushButtonAnimation->setText("Pause simulation");
      _stopwatch.start();
      QMetaObject::invokeMethod(&sim, "step", Qt::QueuedConnection);
      _ui.spinBoxSeed->setEnabled(false);
      _ui.comboBoxDiffusion->setEnabled(false);
      break;
    case SIM_PAUSED:
      if(_simStatus == SIM_RUNNING) sim.lock();
      _pGLWindow->setWindowTitle("Visualization (paused)");
      _ui.pushButtonAnimation->setText("Continue simulation");
      _ui.spinBoxSeed->setEnabled(false);
      _ui.comboBoxDiffusion->setEnabled(false);
      break;
    case SIM_STOPPED:
      if(_simStatus == SIM_RUNNING) sim.lock();
      _pGLWindow->setWindowTitle("Visualization (stopped)");
      _ui.pushButtonAnimation->setText("Continue simulation");
      _ui.spinBoxSeed->setEnabled(true);
      _ui.comboBoxDiffusion->setEnabled(true);
      break;
  }
  _simStatus = newStatus;
}

void MainWindow::on_agentHistButton_clicked(bool checked)
{
  if(checked)
    {
      _pAgentHistogram->updatePlot();
      _pAgentHistogram->show();
    }
  else
    _pAgentHistogram->hide();
}

void MainWindow::saveSettings() const
{
  QSettings settings(QApplication::applicationDirPath() + "/.config", QSettings::IniFormat);
  settings.beginGroup("Simulation");
  settings.setValue("delay", _ui.spinBoxDelay->value());
  settings.setValue("stopClearance", _ui.checkBoxStopClearance->isChecked());
  settings.setValue("stopTime", _ui.checkBoxStopDays->isChecked());
  settings.setValue("stopTimeStep", _ui.spinBoxStopTime->value());
  settings.endGroup();
  settings.beginGroup("Visualization");
  settings.setValue("blend", _ui.checkBoxBlend->isChecked());
  settings.setValue("drawAgents", _ui.checkBoxDrawAgents->isChecked());
  settings.setValue("drawGlyphs", _ui.checkBoxDrawGlyphs->isChecked());
  settings.setValue("drawHeightPlot", _ui.checkBoxDrawHeightPlot->isChecked());
  settings.setValue("drawIsolines", _ui.checkBoxDrawIsolines->isChecked());
  settings.setValue("drawOutcomes", _ui.checkBoxDrawOutcomes->isChecked());
  settings.setValue("drawSmoke", _ui.checkBoxDrawSmoke->isChecked());
  settings.setValue("drawTime", _ui.checkBoxDrawTime->isChecked());
  settings.beginGroup("Granuloma");
  settings.setValue("drawGranulomaBorder", _ui.checkBoxDrawGranulomaBorder->isChecked());
  settings.setValue("granulomaDataset", _ui.comboBoxGranulomaDataset->currentIndex());
  settings.setValue("threshold", _ui.horizontalSliderGranulomaBorderThreshold->value());
  settings.endGroup();
  settings.endGroup();
  settings.beginGroup("Agents");
  _pAgentsWidget->saveSettings();
  settings.endGroup();
  settings.beginGroup("Scalars");
  settings.setValue("smokeDataset", _ui.comboBoxSmokeDataset->currentIndex());
  settings.endGroup();
}
void MainWindow::loadSettings()
{
  QSettings settings(QApplication::applicationDirPath() + "/.config", QSettings::IniFormat);
  settings.beginGroup("Simulation");
  _ui.spinBoxDelay->setValue(settings.value("delay", _ui.spinBoxDelay->value()).toInt());
  _ui.checkBoxStopClearance->setChecked(settings.value("stopClearance", _ui.checkBoxStopClearance->isChecked()).toBool());
  _ui.checkBoxStopDays->setChecked(settings.value("stopTime", _ui.checkBoxStopDays->isChecked()).toBool());
  _ui.spinBoxStopTime->setValue(settings.value("stopTimeStep", _ui.spinBoxStopTime->value()).toBool());
  settings.endGroup();
  settings.beginGroup("Visualization");
  _ui.checkBoxBlend->setChecked(settings.value("blend", _ui.checkBoxBlend->isChecked()).toBool());
  _ui.checkBoxDrawAgents->setChecked(settings.value("drawAgents", _ui.checkBoxDrawAgents->isChecked()).toBool());
  _ui.checkBoxDrawGlyphs->setChecked(settings.value("drawGlyphs", _ui.checkBoxDrawGlyphs->isChecked()).toBool());
  _ui.checkBoxDrawHeightPlot->setChecked(settings.value("drawHeightPlot", _ui.checkBoxDrawHeightPlot->isChecked()).toBool());
  _ui.checkBoxDrawIsolines->setChecked(settings.value("drawIsolines", _ui.checkBoxDrawIsolines->isChecked()).toBool());
  _ui.checkBoxDrawOutcomes->setChecked(settings.value("drawOutcomes", _ui.checkBoxDrawOutcomes->isChecked()).toBool());
  _ui.checkBoxDrawSmoke->setChecked(settings.value("drawSmoke", _ui.checkBoxDrawSmoke->isChecked()).toBool());
  _ui.checkBoxDrawTime->setChecked(settings.value("drawTime", _ui.checkBoxDrawTime->isChecked()).toBool());
  settings.beginGroup("Granuloma");
  _ui.checkBoxDrawGranulomaBorder->setChecked(settings.value("drawGranulomaBorder", _ui.checkBoxDrawGranulomaBorder->isChecked()).toBool());
  _ui.comboBoxGranulomaDataset->setCurrentIndex(settings.value("granulomaDataset", _ui.comboBoxGranulomaDataset->currentIndex()).toInt());
  _ui.horizontalSliderGranulomaBorderThreshold->setValue(settings.value("threshold", _ui.horizontalSliderGranulomaBorderThreshold->value()).toDouble());
  settings.endGroup();
  settings.endGroup();
  _pAgentsWidget->loadSettings();
  settings.beginGroup("Scalars");
  _ui.comboBoxSmokeDataset->setCurrentIndex(settings.value("smokeDataset", _ui.comboBoxSmokeDataset->currentIndex()).toInt());
  settings.endGroup();
}

void MainWindow::setSnapshot(Snapshot* pSnapshot)
{
  delete _pSnapshot;
  _pSnapshot = pSnapshot;
}

void MainWindow::on_groupBoxOutput_toggled(bool arg1)
{
  if(!arg1) {
    if(_pSnapshot) {  //Turn off output
      delete _pSnapshot;
      _pSnapshot = NULL;
    }
    return;
  }
  if (!_pSnapshot)
    {
      QString dirName = QFileDialog::getExistingDirectory(this, "Browse directory");
      if (dirName == QString::null) {
        _ui.groupBoxOutput->setChecked(false);
        return;
      }
      _ui.lineEditOutput->setText(dirName);
      QString fileName = QString("%1%2seed%3.csv").arg(dirName).arg(QDir::separator()).arg(g_Rand.getSeed());
      const char* strFileName = fileName.toLatin1().data();

      _pSnapshot = new Snapshot(dirName, strFileName);

      if (!_pSnapshot->isGood())
        {
          QMessageBox::critical(this, "Lung ABM",
                                "Failed to open/create file '" + fileName + "' for writing.",
                                QMessageBox::Ok, QMessageBox::Ok);
          delete _pSnapshot;
          _pSnapshot = NULL;
          _ui.groupBoxOutput->setChecked(false);
        }
    }
}

void MainWindow::on_comboBoxColorMapSource_currentIndexChanged(int index)
{

  refreshCurrentColorMap(COLORMAP_SRC(index));
  refreshCurrentNormalizer(COLORMAP_SRC(index));
  refreshCurrentScalarGrid(COLORMAP_SRC(index));
  Q_CHECK_PTR(_pCurrentColorMap);
  Q_CHECK_PTR(_pCurrentNormalizer);

  float min = _pCurrentNormalizer->getMin();
  float max = _pCurrentNormalizer->getMax();

  if (_pCurrentNormalizer->getClamping())
      _ui.comboBoxMappingMethod->setCurrentIndex(0);
  else
      _ui.comboBoxMappingMethod->setCurrentIndex(1);

  _ui.doubleSpinBoxMin->setValue(min);
  _ui.doubleSpinBoxMax->setValue(max);

  int hue = (int)(_pCurrentColorMap->getHueDelta() * _SLIDER_COLORMAP_HUE);
  int sat = (int)(_pCurrentColorMap->getSatDelta() * _SLIDER_COLORMAP_SAT);
  int val = (int)(_pCurrentColorMap->getValDelta() * _SLIDER_COLORMAP_VAL);
  int alpha = (int)(_pCurrentColorMap->getAlpha() * _SLIDER_COLORMAP_ALPHA);

  bool reverse = _pCurrentColorMap->getInvert();
  int nrBands = _pCurrentColorMap->getNrBands();

  _ui.horizontalSliderHue->setValue(hue);
  _ui.horizontalSliderSat->setValue(sat);
  _ui.horizontalSliderVal->setValue(val);
  _ui.horizontalSliderAlpha->setValue(alpha);

  COLORMAP idx = _pCurrentColorMap->getType();
  _ui.comboBoxScalarColoring->setCurrentIndex(idx);

  _ui.checkBoxReverseColorMap->setChecked(reverse);
  _ui.spinBoxNumberOfColors->setValue(nrBands);

  setColorMap(idx);
}

void MainWindow::on_comboBoxSmokeDataset_currentIndexChanged(int index)
{
  _pItfc->setScalarSmokeDataset(getNewScalarDataset(index));

  if (_ui.comboBoxMappingMethod->currentIndex() == 1 &&
      _ui.comboBoxColorMapSource->currentIndex() == SMOKE_SRC)
    {
      rescaleRequest();
    }
  _pItfc->updateGrids();

  emit updateGL();
}

void MainWindow::on_comboBoxGlyphScalarDataset_currentIndexChanged(int index)
{
  _pItfc->setScalarGlyphDataset(getNewScalarDataset(index));

  if (_ui.comboBoxMappingMethod->currentIndex() == 1 &&
      _ui.comboBoxColorMapSource->currentIndex() == GLYPHS_SRC)
    {
      rescaleRequest();
    }
  _pItfc->updateGrids();

  emit updateGL();
}

void MainWindow::on_comboBoxGlyphVectorDataset_currentIndexChanged(int index)
{
  _pItfc->setScalarGlyphDataset(getNewScalarDataset(index));

  if (_ui.comboBoxMappingMethod->currentIndex() == 1 &&
      _ui.comboBoxColorMapSource->currentIndex() == GLYPHS_SRC)
    {
      rescaleRequest();
    }
  _pItfc->updateGrids();

  emit updateGL();
}

void MainWindow::on_doubleSpinBoxGlyphScaling_valueChanged(double arg1)
{
  _pItfc->setGlyphScalingFactor(arg1);
  emit updateGL();
}

void MainWindow::on_checkBoxDrawAgents_toggled(bool checked)
{
  _pItfc->setDrawAgents(checked);
  _ui.tabAgents->setEnabled(checked);

  emit updateGL();
}

void MainWindow::on_checkBoxDrawSmoke_toggled(bool checked)
{
  _pItfc->setDrawSmoke(checked);
  _ui.tabSmoke->setEnabled(checked || _pItfc->getDrawIsolines());

  _pItfc->updateGrids();
  emit updateGL();
}

void MainWindow::on_checkBoxDrawGlyphs_toggled(bool checked)
{
  _pItfc->setDrawGlyphs(checked);
  _ui.tabGlyphs->setEnabled(checked);

  _pItfc->updateGrids();
  emit updateGL();
}

void MainWindow::on_checkBoxDrawIsolines_toggled(bool checked)
{
  _pItfc->setDrawIsolines(checked);
  _ui.groupBoxIsoline->setEnabled(checked);
  _ui.tabSmoke->setEnabled(checked || _pItfc->getDrawSmoke());

  _pItfc->updateGrids();
  emit updateGL();
}

void MainWindow::on_checkBoxDrawHeightPlot_toggled(bool checked)
{
  if (checked)
    {
      _ui.checkBoxDrawGlyphs->setChecked(false);
      _ui.checkBoxDrawSmoke->setChecked(false);
      _ui.checkBoxDrawIsolines->setChecked(false);
      _ui.checkBoxDrawGlyphs->setEnabled(false);
      _ui.checkBoxDrawSmoke->setEnabled(false);
      _ui.checkBoxDrawIsolines->setEnabled(false);
    }
  else
    {
      _ui.checkBoxDrawGlyphs->setEnabled(true);
      _ui.checkBoxDrawSmoke->setEnabled(true);
      _ui.checkBoxDrawIsolines->setEnabled(true);
    }

  _pItfc->setDrawHeightPlot(checked);
  _ui.tabHeightPlot->setEnabled(checked);

  _pItfc->updateGrids();
  emit updateGL();
}

void MainWindow::on_checkBoxBlend_toggled(bool checked)
{
  _pItfc->setBlend(checked);
  _ui.labelTransparency->setEnabled(checked);
  _ui.horizontalSliderAlpha->setEnabled(checked);
  _ui.labelHeightGridTransparency->setEnabled(checked);
  _ui.horizontalSliderHeightAlpha->setEnabled(checked);

  emit updateGL();
}

void MainWindow::on_checkBoxGlyphClamping_toggled(bool checked)
{
  _pItfc->setGlyphClamping(checked);
  emit updateGL();
}
