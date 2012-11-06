#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QTime>
#include <QtGui/QMainWindow>
#include <QMessageBox>
#include <QString>
#include <QDir>
#include "ui_mainwindow.h"

#include "gui/agenthistogram.h"
#include "maininterface.h"


class GLWindow;
class StatWidget;
class AgentsWidget;
class GraphController;
class Snapshot;

enum SimStatus { SIM_RUNNING, SIM_PAUSED, SIM_STOPPED };


enum COLORMAP_SRC { SMOKE_SRC=0, GLYPHS_SRC, ISOLINES_SRC, HEIGHTPLOT_SRC, NCOLORMAPSRCS };

inline std::ostream& operator<<(std::ostream& s, COLORMAP_SRC m) {
  switch(m) {
    case SMOKE_SRC:      s<<"Smoke";       break;
    case GLYPHS_SRC:     s<<"Glyphs";      break;
    case ISOLINES_SRC:   s<<"Isolines";    break;
    case HEIGHTPLOT_SRC: s<<"Height plot"; break;
    default: assert(!"Invalid colormap source"); break;
  }
  return s;
}

class MainWindow : public QMainWindow
{
  Q_OBJECT

public:
  MainWindow(MainInterface* pItfc, GLWindow* pGLWindow, QWidget* pParamWindow,
             StatWidget* pStatWidget, AgentsWidget* pAgentsWidget, const QDir& dir, QWidget* parent = 0);
  ~MainWindow();
  Ui::MainWindowClass& getUI();
  void setScriptingMode(bool scriptingMode);
  void setSnapshot(Snapshot* pSnapshot);
  void loadState(std::string fileName);

  static int getScalarDataSetIndex(const QString &dataSetName);
  static int getScalarDataSetIndex(const std::string &dataSetName);
  static void getScalarDataSetNames(std::string &names);

  static MainWindow* _pMainWindow;
  static const QString _COLORMAP_GRAYSCALE;
  static const QString _COLORMAP_RAINBOW;
  static const QString _COLORMAP_FIRE;
  static const QString _COLORMAP_COOLWARM;
  static const QString _COLORMAP_GREENRED;

  static const QString _MAP_METHOD_SCALE;
  static const QString _MAP_METHOD_CLAMP;

  static const QString _DATASET_TNF;
  static const QString _DATASET_TNF_GRADIENT;
  static const QString _DATASET_CCL2;
  static const QString _DATASET_CCL5_GRADIENT;
  static const QString _DATASET_CCL5;
  static const QString _DATASET_CCL2_GRADIENT;
  static const QString _DATASET_CXCL9;
  static const QString _DATASET_CXCL9_GRADIENT;
  static const QString _DATASET_IL10;
  static const QString _DATASET_IL10_GRADIENT;
  static const QString _DATASET_ATTRACTANT;
  static const QString _DATASET_ATTRACTANT_GRADIENT;
  static const QString _DATASET_TNF_ATTR_EXTMTB;
  static const QString _DATASET_CELL_DENSITY;
  static const QString _DATASET_EXTMTB;
  static const QString _DATASET_INTMTB;
  static const QString _DATASET_TOTMTB;
  static const QString _DATASET_KILLINGS;
  static const QString _DATASET_KILLINGS_GRADIENT;

  static const QString _SCALAR_DATASETS[];
  static const int _N_SCALAR_DATASETS;


  static const QString _INTERPOLATION_NEAREST_NEIGHBOR;
  static const QString _INTERPOLATION_BILINEAR;

  static const QString _COLORMAP_SOURCE_SMOKE;
  static const QString _COLORMAP_SOURCE_GLYPHS;
  static const QString _COLORMAP_SOURCE_ISOLINES;
  static const QString _COLORMAP_SOURCE_HEIGHTPLOT;

  static const QString _GLYPH_HEDGEHOG;
  static const QString _GLYPH_CONE;
  static const QString _GLYPH_ARROW;
  static const QString _GLYPH_TEXTURE;

  static const QString _DIFF_SOR_WRONG;
  static const QString _DIFF_SOR_CORRECT;
  static const QString _DIFF_FTCS;
  static const QString _DIFF_FTCS_SWAP;
  static const QString _DIFF_ADE_SWAP;

  static const QString _OUTCOME_MTB;
  static const QString _OUTCOME_AREA;
  static const QString _OUTCOME_NONE;
  static const QString _OUTCOME_METHOD_1;
  static const QString _OUTCOME_METHOD_2;

  static const int _TIMER_PERIOD = 15; // In milliseconds (60 Hz == 16.6ms)

public slots:
  void toggleAnimation();
  void setBlend();
  void setDrawAgents(bool checked);
  void setDrawSmoke(bool checked);
  void setDrawGlyphs(bool checked);
  void setDrawIsolines(bool checked);
  void setDrawHeightPlot(bool checked);
  void setGlyphScaling(int value);
  void setColorMap(int value);
  void setMappingMethod(int value);
  void setSmokeDataset(const QString& value);
  void updateColorMap();
  void updateColorMapLabels();
  void updateHeightMinMax();
  void rescaleRequest();
  void rescaleHeightRequest();
  void updateIsolinesSettings();
  void updateHeightPlotSettings();
  void setGlyphVectorDataset(const QString& value);
  void setGlyphScalarDataset(const QString& value);
  void setGlyph(const QString& value);
  void setGlyphClamping(bool checked);
  void dumpGrids();
  void selectSmokeColorMap();
  void selectGlyphsColorMap();
  void selectIsolinesColorMap();
  void selectHeightPlotColorMap();
  void selectColorMap0();
  void selectColorMap1();
  void selectColorMap2();
  void selectColorMap3();
  void selectColorMap4();
  void updateGranulomaSettings(void);
  void updateOutcomeSettings();
  void updateOutcomeParameters();
  void stop();
  void loadState();
  void saveState();
  void takePictureSnapshot();
  void saveSettings() const;
  void loadSettings();

signals:
  void updateGL();
  void updateSelectedCellStats();
  void updateColorMap(ColorMap* pColorMap);
  void updateColorMapLabels(float min, float max);

protected:
  void timerEvent(QTimerEvent* pEvent);
  void closeEvent(QCloseEvent* pEvent);

private slots:
  void on_agentHistButton_clicked(bool checked);
  void update();

  void on_checkBoxStopDays_toggled(bool checked);

  void on_spinBoxStopTime_valueChanged(int arg1);

  void on_spinBoxStopDays_valueChanged(int arg1);

  void on_checkBoxStopClearance_toggled(bool checked);

  void on_pushButtonShowParams_clicked();

  void on_actionShowParameters_triggered();

  void on_comboBoxDiffusion_currentIndexChanged(int index);

  void on_spinBoxDelay_valueChanged(int arg1);

  void on_spinBoxSeed_valueChanged(int arg1);

  void on_groupBoxOutput_toggled(bool arg1);

  void on_comboBoxColorMapSource_currentIndexChanged(int index);

private:
  void initVisualizationTab();
  void initSimulationTab();
  void initOutcomeTab();
  void initColorMapTab();
  void initSmokeTab();
  void initGlyphsTab();
  void initIsolinesTab();
  void initHeightTab();
  void initOutputTab();
  void initShortcuts();
  void initComboAllScalars(QComboBox *&comboBox);
  void refreshCurrentNormalizer(const COLORMAP_SRC m);
  void refreshCurrentColorMap(const COLORMAP_SRC colorMapSource);
  void refreshCurrentScalarGrid(const COLORMAP_SRC m);
  ScalarDataset* getNewScalarDataset(const QString& value);
  VectorDataset* getNewVectorDataset(const QString& value);
  void updateStatLabels(const Stats& stats);
  void switchStatus(SimStatus newStatus);
  void updateWindowStatus();
  bool isTimeForAction(int actionInterval, int simTime);

  Ui::MainWindowClass _ui;
  MainInterface* _pItfc;
  AgentHistogram* _pAgentHistogram;
  GraphController* _pGraphController;
  GLWindow* _pGLWindow;
  QWidget* _pParamWindow;
  StatWidget* _pStatWidget;
  AgentsWidget* _pAgentsWidget;
  ColorMap* _pCurrentColorMap;
  ScalarNormalizer* _pCurrentNormalizer;
  ScalarGrid* _pCurrentScalarGrid;
  int _timerId;
  SimStatus _simStatus;
  QTime _stopwatch;
  Snapshot* _pSnapshot;
  bool _scriptingMode;
  QDir _dir;
};

inline void MainWindow::selectSmokeColorMap()
{
  _ui.comboBoxColorMapSource->setCurrentIndex(0);
}

inline void MainWindow::selectGlyphsColorMap()
{
  _ui.comboBoxColorMapSource->setCurrentIndex(1);
}

inline void MainWindow::selectIsolinesColorMap()
{
  _ui.comboBoxColorMapSource->setCurrentIndex(2);
}

inline void MainWindow::selectHeightPlotColorMap()
{
  _ui.comboBoxColorMapSource->setCurrentIndex(3);
}

inline void MainWindow::selectColorMap0()
{
  _ui.comboBoxScalarColoring->setCurrentIndex(0);
}

inline void MainWindow::selectColorMap1()
{
  _ui.comboBoxScalarColoring->setCurrentIndex(1);
}

inline void MainWindow::selectColorMap2()
{
  _ui.comboBoxScalarColoring->setCurrentIndex(2);
}

inline void MainWindow::selectColorMap3()
{
  _ui.comboBoxScalarColoring->setCurrentIndex(3);
}
inline void MainWindow::selectColorMap4()
{
  _ui.comboBoxScalarColoring->setCurrentIndex(4);
}

inline void MainWindow::refreshCurrentColorMap(const COLORMAP_SRC m)
{
  switch(m) {
    case SMOKE_SRC:      _pCurrentColorMap = _pItfc->getColorMapSmoke(); break;
    case GLYPHS_SRC:     _pCurrentColorMap = _pItfc->getColorMapGlyphs(); break;
    case ISOLINES_SRC:   _pCurrentColorMap = _pItfc->getColorMapIsolines(); break;
    case HEIGHTPLOT_SRC: _pCurrentColorMap = _pItfc->getColorMapHeightPlot(); break;
    default: assert(!"Invalid colormap src"); break;
  }
}

inline void MainWindow::refreshCurrentNormalizer(const COLORMAP_SRC m)
{
  switch(m) {
    case SMOKE_SRC:      _pCurrentNormalizer = _pItfc->getScalarSmokeNormalizer(); break;
    case GLYPHS_SRC:     _pCurrentNormalizer = _pItfc->getScalarGlyphNormalizer(); break;
    case ISOLINES_SRC:   _pCurrentNormalizer = _pItfc->getScalarIsolinesNormalizer(); break;
    case HEIGHTPLOT_SRC: _pCurrentNormalizer = _pItfc->getScalarHeightColorNormalizer(); break;
    default: assert(!"Invalid colormap src"); break;
  }
}

inline void MainWindow::refreshCurrentScalarGrid(const COLORMAP_SRC m)
{
  switch(m) {
    case SMOKE_SRC:
    case ISOLINES_SRC:   _pCurrentScalarGrid = _pItfc->getScalarSmokeGrid(); break;
    case GLYPHS_SRC:     _pCurrentScalarGrid = _pItfc->getScalarGlyphGrid(); break;
    case HEIGHTPLOT_SRC: _pCurrentScalarGrid = _pItfc->getScalarHeightColorGrid(); break;
    default: assert(!"Invalid colormap src"); break;
  }
}

inline Ui::MainWindowClass& MainWindow::getUI()
{
  return _ui;
}

inline void MainWindow::setScriptingMode(bool scriptingMode)
{
  _scriptingMode = scriptingMode;
}

#endif // MAINWINDOW_H
