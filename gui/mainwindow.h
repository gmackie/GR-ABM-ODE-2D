#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QTime>
#include <QtGui/QMainWindow>
#include <QMessageBox>
#include <QString>
#include "ui_mainwindow.h"

#include "gui/glwindow.h"
#include "gui/paramwindow.h"
#include "gui/statwidget.h"
#include "gui/agentswidget.h"
#include "maininterface.h"
#include "snapshot.h"
#include "colormaps/colormap.h"
#include "scalardatasets/scalarnormalizer.h"
#include "scalardatasets/scalargrid.h"

typedef enum { SIM_RUNNING, SIM_PAUSED, SIM_STOPPED } SimStatus;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(MainInterface* pItfc, GLWindow* pGLWindow, QWidget* pParamWindow,
    		StatWidget* pStatWidget, QWidget* pAgentsWidget, QWidget* parent = 0);
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

    static const QString _OUTCOME_MTB;
    static const QString _OUTCOME_AREA;
    static const QString _OUTCOME_NONE;
    static const QString _OUTCOME_METHOD_1;
    static const QString _OUTCOME_METHOD_2;

    static const int _TIMER_PERIOD = 10; // In milliseconds.

public slots:
	void updateTimeBox(int);
	void toggleAnimation();
	void setBlend();
	void setDrawAgents(bool checked);
	void setDrawSmoke(bool checked);
	void setDrawGlyphs(bool checked);
	void setDrawIsolines(bool checked);
	void setDrawHeightPlot(bool checked);
	void setGlyphScaling(int value);
	void setColorMap(const QString& value);
	void setMappingMethod(const QString& value);
	void setSmokeDataset(const QString& value);
	void setColorMapSource(const QString& value);
	void updateDelay(int value);
	void updateSeed(int value);
	void updateColorMap();
	void updateColorMapLabels();
	void updateHeightMinMax();
	void updateDiffusionMethod(const QString& value);
	void updateStopCriteria();
	void rescaleRequest();
	void rescaleHeightRequest();
	void updateIsolinesSettings();
	void updateHeightPlotSettings();
	void setGlyphVectorDataset(const QString& value);
	void setGlyphScalarDataset(const QString& value);
	void setGlyph(const QString& value);
	void setGlyphClamping(bool checked);
	void dumpGrids();
	void showParams();
	void selectSmokeColorMap();
	void selectGlyphsColorMap();
	void selectIsolinesColorMap();
	void selectHeightPlotColorMap();
	void selectColorMap0();
	void selectColorMap1();
	void selectColorMap2();
	void selectColorMap3();
	void selectColorMap4();
	void setEnableOutput(bool checked);
	void updateGranulomaSettings(void);
	void updateOutcomeSettings();
	void updateOutcomeParameters();
	void stop();
	void loadState();
	void saveState();
	void takePictureSnapshot();

signals:
	void updateGL();
	void updateSelectedCellStats();
	void updateColorMap(ColorMap* pColorMap);
	void updateColorMapLabels(float min, float max);

protected:
	void timerEvent(QTimerEvent* pEvent);
	void closeEvent(QCloseEvent* pEvent);

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
	void refreshCurrentNormalizer(const QString& colorMapSource);
	void refreshCurrentColorMap(const QString& colorMapSource);
	void refreshCurrentScalarGrid(const QString& colorMapSource);
	ScalarDataset* getNewScalarDataset(const QString& value);
	VectorDataset* getNewVectorDataset(const QString& value);
	void updateStatLabels(const GrStat& stats);
	void switchStatus(SimStatus newStatus);
	void updateWindowStatus();
	bool isTimeForAction(int actionInterval, int simTime);

	Ui::MainWindowClass _ui;
	MainInterface* _pItfc;
	GLWindow* _pGLWindow;
	QWidget* _pParamWindow;
	StatWidget* _pStatWidget;
	QWidget* _pAgentsWidget;
	ColorMap* _pCurrentColorMap;
    ScalarNormalizer* _pCurrentNormalizer;
    ScalarGrid* _pCurrentScalarGrid;
	int _timerId;
	SimStatus _simStatus;
	QTime _stopwatch;
	Snapshot* _pSnapshot;
	bool _scriptingMode;
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

inline void MainWindow::refreshCurrentColorMap(const QString& colorMapSource)
{
	if (colorMapSource == _COLORMAP_SOURCE_SMOKE)
	{
		_pCurrentColorMap = _pItfc->getColorMapSmoke();
	}
	else if (colorMapSource == _COLORMAP_SOURCE_GLYPHS)
	{
		_pCurrentColorMap = _pItfc->getColorMapGlyphs();
	}
	else if (colorMapSource == _COLORMAP_SOURCE_ISOLINES)
	{
		_pCurrentColorMap = _pItfc->getColorMapIsolines();
	}
	else if (colorMapSource == _COLORMAP_SOURCE_HEIGHTPLOT)
	{
		_pCurrentColorMap = _pItfc->getColorMapHeightPlot();
	}
	else
	{
		assert(false);
	}
}

inline void MainWindow::refreshCurrentNormalizer(const QString& colorMapSource)
{
	if (colorMapSource == _COLORMAP_SOURCE_SMOKE)
	{
		_pCurrentNormalizer = _pItfc->getScalarSmokeNormalizer();
	}
	else if (colorMapSource == _COLORMAP_SOURCE_GLYPHS)
	{
		_pCurrentNormalizer = _pItfc->getScalarGlyphNormalizer();
	}
	else if (colorMapSource == _COLORMAP_SOURCE_ISOLINES)
	{
		_pCurrentNormalizer = _pItfc->getScalarIsolinesNormalizer();
	}
	else if (colorMapSource == _COLORMAP_SOURCE_HEIGHTPLOT)
	{
		_pCurrentNormalizer = _pItfc->getScalarHeightColorNormalizer();
	}
	else
	{
		assert(false);
	}
}

inline void MainWindow::refreshCurrentScalarGrid(const QString& colorMapSource)
{
	if (colorMapSource == _COLORMAP_SOURCE_SMOKE ||
		colorMapSource == _COLORMAP_SOURCE_ISOLINES)
	{
		_pCurrentScalarGrid = _pItfc->getScalarSmokeGrid();
	}
	else if (colorMapSource == _COLORMAP_SOURCE_GLYPHS)
	{
		_pCurrentScalarGrid = _pItfc->getScalarGlyphGrid();
	}
	else if (colorMapSource == _COLORMAP_SOURCE_HEIGHTPLOT)
	{
		_pCurrentScalarGrid = _pItfc->getScalarHeightColorGrid();
	}
	else
	{
		assert(false);
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

inline void MainWindow::setSnapshot(Snapshot* pSnapshot)
{
	delete _pSnapshot;
	_pSnapshot = pSnapshot;
}

#endif // MAINWINDOW_H
