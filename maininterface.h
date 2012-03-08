/*
 * maininterface.h
 *
 *  Created on: 21-sep-2008
 *      Author: S030858
 */

#ifndef MAININTERFACE_H_
#define MAININTERFACE_H_

#include "scalardatasets/scalaragentgridbase.h"
#include "scalardatasets/scalarnormalizer.h"
#include "scalardatasets/scalardataset.h"
#include "scalardatasets/scalargrid.h"
#include "vectordatasets/vectordataset.h"
#include "vectordatasets/vectorgrid.h"
#include "vectordatasets/vectorgradientdataset.h"
#include "vectordatasets/vector.h"
#include "simulation.h"
#include "visualization/agentsvisualization.h"
#include "visualization/invisiblequadvisualisation.h"
#include "visualization/smokevisualization.h"
#include "visualization/glyphvisualization.h"
#include "visualization/isolinesvisualization.h"
#include "visualization/heightplotvisualization.h"
#include "glyphs/glyph.h"
#include "glyphs/glyphhedgehog.h"
#include "grviz.h"
#include <fstream>

class MainInterface
{
public:
	MainInterface(const Pos& dim, Visualization* pAgentVisualization, ScalarAgentGridBase* pScalarAgentGrid);
	virtual ~MainInterface();
	void doStep();
	ScalarNormalizer* getScalarSmokeNormalizer();
	ScalarNormalizer* getScalarGlyphNormalizer();
	ScalarNormalizer* getScalarIsolinesNormalizer();
	ScalarNormalizer* getScalarHeightNormalizer();
	ScalarNormalizer* getScalarHeightColorNormalizer();
	ScalarDataset* getScalarSmokeDataset();
	void setScalarSmokeDataset(ScalarDataset* pScalarDataset);
	ScalarDataset* getScalarGranulomaDataset();
	void setScalarGranulomaDataset(ScalarDataset* pScalarDataset);
	ScalarDataset* getScalarGlyphDataset();
	void setScalarGlyphDataset(ScalarDataset* pScalarDataset);
	VectorDataset* getVectorGlyphDataset();
	void setVectorGlyphDataset(VectorDataset* pVectorDataset);
	VectorDataset* getVectorSmokeDataset();
	void setVectorSmokeDataset(VectorDataset* pVectorDataset);
	ScalarDataset* getScalarHeightDataset();
	void setScalarHeightDataset(ScalarDataset* pScalarDataset);
	ScalarDataset* getScalarHeightColorDataset();
	void setScalarHeightColorDataset(ScalarDataset* pScalarDataset);
	ColorMap* getColorMapSmoke();
	void setColorMapSmoke(ColorMap* pColorMap);
	ColorMap* getColorMapGlyphs();
	void setColorMapGlyphs(ColorMap* pColorMap);
	ColorMap* getColorMapIsolines();
	void setColorMapIsolines(ColorMap* pColorMap);
	ColorMap* getColorMapHeightPlot();
	void setColorMapHeightPlot(ColorMap* pColorMap);
	bool getDrawAgents();
	void setDrawAgents(bool drawAgents);
	bool getDrawSmoke();
	void setDrawSmoke(bool drawSmoke);
	bool getDrawGlyphs();
	void setDrawGlyphs(bool drawGlyphs);
	bool getDrawIsolines();
	void setDrawIsolines(bool drawIsolines);
	bool getDrawHeightPlot();
	void setDrawHeightPlot(bool drawHeightPlot);
	bool getDrawGranuloma();
	void setDrawGranuloma(bool drawGranuloma);
	Simulation& getSimulation();
	void visualize();
	ScalarAgentGridBase* getScalarAgentGrid();
	ScalarGrid* getScalarSmokeGrid();
	ScalarGrid* getScalarGlyphGrid();
	VectorGrid* getVectorGlyphGrid();
	ScalarGrid* getScalarHeightGrid();
	ScalarGrid* getScalarHeightColorGrid();
	void setGlyphGridNN(bool value);
	void setSmokeGridNN(bool value);
	void setHeightGridNN(bool value);
	void setGlyphScalingFactor(int scaleFactor);
	void updateGrids();
	void setGlyph(Glyph* pGlyph);
	void setGlyphClamping(bool clamp);
	const std::vector<float>& getIsolinesValueVector();
	void setIsolinesValueVector(const std::vector<float>& vec);
	float getIsolinesWidth();
	void setIsolinesWidth(float width);
	float getHeightMax();
	void setHeightMax(float height);
	float getGridHeightMax();
	void setGridHeightMax(float height);
	bool getDrawHeightGrid();
	void setDrawHeightGrid(bool value);
	void serializeGrids(std::ofstream& outFile);
    bool getBlend();
    void setBlend(bool blend);
    void setHeightGridAlpha(float alpha);
    void setGranulomaThreshold(float refValue);
    int getTime();
    void setTime(int time);
    void incTime(int dTime);
    int getSimTime();
    void setSimTime(int simTime);

private:
    Simulation* _simulation;
    ScalarNormalizer _scalarSmokeNormalizer;
    ScalarNormalizer _scalarGlyphNormalizer;
    ScalarNormalizer _scalarIsolinesNormalizer;
    ScalarNormalizer _scalarHeightNormalizer;
    ScalarNormalizer _scalarHeightColorNormalizer;
    ScalarNormalizer _scalarGranulomaNormalizer;
    ScalarDataset* _pScalarSmokeDataset;
    ScalarDataset* _pScalarGlyphDataset;
    VectorDataset* _pVectorGlyphDataset;
    ScalarDataset* _pScalarHeightDataset;
    ScalarDataset* _pScalarHeightColorDataset;
    ScalarDataset* _pScalarGranulomaDataset;
    ScalarAgentGridBase* _pScalarAgentGrid;
    ScalarGrid _scalarSmokeGrid;
    ScalarGrid _scalarGlyphGrid;
    VectorGrid _vectorGlyphGrid;
    ScalarGrid _scalarHeightGrid;
    ScalarGrid _scalarHeightColorGrid;
    ScalarGrid _scalarGranulomaGrid;
    Visualization* _pAgentsVisualization;
    InvisibleQuadVisualisation _invisibleQuadVisualization;
    SmokeVisualization _smokeVisualization;
    GlyphVisualization _glyphVisualization;
    IsolinesVisualization _isolinesVisualization;
    HeightPlotVisualization _heightPlotVisualization;
    IsolinesVisualization _granulomaVisualization;
    ColorMap* _pColorMapSmoke;
    ColorMap* _pColorMapGlyphs;
    ColorMap* _pColorMapIsolines;
    ColorMap* _pColorMapHeightPlot;
    ColorMap* _pColorMapGranuloma;
    bool _drawAgents;
    bool _drawSmoke;
    bool _drawGlyphs;
    bool _drawIsolines;
    bool _drawHeightPlot;
    bool _drawGranuloma;
    bool _glyphGridNN;
    bool _smokeGridNN;
    bool _heightGridNN;
    bool _blend;
	int _time;
	int _simTime;
};

inline void MainInterface::doStep()
{
	updateGrids();
}

inline ScalarNormalizer* MainInterface::getScalarSmokeNormalizer()
{
	return &_scalarSmokeNormalizer;
}

inline ScalarNormalizer* MainInterface::getScalarGlyphNormalizer()
{
	return &_scalarGlyphNormalizer;
}

inline ScalarNormalizer* MainInterface::getScalarIsolinesNormalizer()
{
	return &_scalarIsolinesNormalizer;
}

inline ScalarNormalizer* MainInterface::getScalarHeightNormalizer()
{
	return &_scalarHeightNormalizer;
}

inline ScalarNormalizer* MainInterface::getScalarHeightColorNormalizer()
{
	return &_scalarHeightColorNormalizer;
}

inline ScalarDataset* MainInterface::getScalarSmokeDataset()
{
	return _pScalarSmokeDataset;
}

inline void MainInterface::setScalarSmokeDataset(ScalarDataset* pScalarDataset)
{
	delete _pScalarSmokeDataset;
	_pScalarSmokeDataset = pScalarDataset;
}

inline ScalarDataset* MainInterface::getScalarGranulomaDataset()
{
	return _pScalarGranulomaDataset;
}

inline void MainInterface::setScalarGranulomaDataset(ScalarDataset* pScalarDataset)
{
	delete _pScalarGranulomaDataset;
	_pScalarGranulomaDataset = pScalarDataset;
}


inline ScalarDataset* MainInterface::getScalarGlyphDataset()
{
	return _pScalarGlyphDataset;
}

inline void MainInterface::setScalarGlyphDataset(ScalarDataset* pScalarDataset)
{
	delete _pScalarGlyphDataset;
	_pScalarGlyphDataset = pScalarDataset;
}

inline VectorDataset* MainInterface::getVectorGlyphDataset()
{
	return _pVectorGlyphDataset;
}

inline void MainInterface::setVectorGlyphDataset(VectorDataset* pVectorDataset)
{
	delete _pVectorGlyphDataset;
	_pVectorGlyphDataset = pVectorDataset;
}

inline ScalarDataset* MainInterface::getScalarHeightDataset()
{
	return _pScalarHeightDataset;
}

inline void MainInterface::setScalarHeightDataset(ScalarDataset* pScalarDataset)
{
	delete _pScalarHeightDataset;
	_pScalarHeightDataset = pScalarDataset;

	VectorGradientDataset* pGradientDataset =
		new VectorGradientDataset(_pScalarHeightDataset, true);

	_heightPlotVisualization.setGradientDataset(pGradientDataset);
}

inline ScalarDataset* MainInterface::getScalarHeightColorDataset()
{
	return _pScalarHeightColorDataset;
}

inline void MainInterface::setScalarHeightColorDataset(ScalarDataset* pScalarDataset)
{
	delete _pScalarHeightColorDataset;
	_pScalarHeightColorDataset = pScalarDataset;
}

inline ColorMap* MainInterface::getColorMapSmoke()
{
	return _pColorMapSmoke;
}

inline void MainInterface::setColorMapSmoke(ColorMap* pColorMap)
{
	delete _pColorMapSmoke;
	_pColorMapSmoke = pColorMap;
}

inline ColorMap* MainInterface::getColorMapGlyphs()
{
	return _pColorMapGlyphs;
}

inline void MainInterface::setColorMapGlyphs(ColorMap* pColorMap)
{
	delete _pColorMapGlyphs;
	_pColorMapGlyphs = pColorMap;
}

inline ColorMap* MainInterface::getColorMapIsolines()
{
	return _pColorMapIsolines;
}

inline void MainInterface::setColorMapIsolines(ColorMap* pColorMap)
{
	delete _pColorMapIsolines;
	_pColorMapIsolines = pColorMap;
}

inline ColorMap* MainInterface::getColorMapHeightPlot()
{
	return _pColorMapHeightPlot;
}

inline void MainInterface::setColorMapHeightPlot(ColorMap* pColorMap)
{
	delete _pColorMapHeightPlot;
	_pColorMapHeightPlot = pColorMap;
}

inline bool MainInterface::getDrawAgents()
{
	return _drawAgents;
}

inline void MainInterface::setDrawAgents(bool drawAgents)
{
	_drawAgents = drawAgents;
}

inline bool MainInterface::getDrawSmoke()
{
	return _drawSmoke;
}

inline void MainInterface::setDrawSmoke(bool drawSmoke)
{
	_drawSmoke = drawSmoke;
}

inline bool MainInterface::getDrawGlyphs()
{
	return _drawGlyphs;
}

inline void MainInterface::setDrawGlyphs(bool drawGlyphs)
{
	_drawGlyphs = drawGlyphs;
}

inline bool MainInterface::getDrawIsolines()
{
	return _drawIsolines;
}

inline void MainInterface::setDrawIsolines(bool drawIsolines)
{
	_drawIsolines = drawIsolines;
}

inline bool MainInterface::getDrawHeightPlot()
{
	return _drawHeightPlot;
}

inline void MainInterface::setDrawHeightPlot(bool drawHeightPlot)
{
	_drawHeightPlot = drawHeightPlot;
}

inline bool MainInterface::getDrawGranuloma()
{
	return _drawGranuloma;
}

inline void MainInterface::setDrawGranuloma(bool drawGranuloma)
{
	_drawGranuloma = drawGranuloma;
}

inline Simulation& MainInterface::getSimulation()
{
	return *_simulation;
}

inline ScalarAgentGridBase* MainInterface::getScalarAgentGrid()
{
	return _pScalarAgentGrid;
}

inline ScalarGrid* MainInterface::getScalarSmokeGrid()
{
	return &_scalarSmokeGrid;
}

inline ScalarGrid* MainInterface::getScalarGlyphGrid()
{
	return &_scalarGlyphGrid;
}

inline VectorGrid* MainInterface::getVectorGlyphGrid()
{
	return &_vectorGlyphGrid;
}

inline ScalarGrid* MainInterface::getScalarHeightGrid()
{
	return &_scalarHeightGrid;
}

inline ScalarGrid* MainInterface::getScalarHeightColorGrid()
{
	return &_scalarHeightColorGrid;
}

inline void MainInterface::setGlyphGridNN(bool value)
{
	_glyphGridNN = value;
}

inline void MainInterface::setSmokeGridNN(bool value)
{
	_smokeGridNN = value;
}

inline void MainInterface::setHeightGridNN(bool value)
{
	_heightGridNN = value;
}

inline void MainInterface::setGlyphScalingFactor(int scaleFactor)
{
	_glyphVisualization.setScaleFactor(scaleFactor);
}

inline void MainInterface::setGlyph(Glyph* pGlyph)
{
	_glyphVisualization.setGlyph(pGlyph);
}

inline void MainInterface::setGlyphClamping(bool clamp)
{
	_glyphVisualization.setClamping(clamp);
}

inline const std::vector<float>& MainInterface::getIsolinesValueVector()
{
	return _isolinesVisualization.getTargetValueVector();
}

inline void MainInterface::setIsolinesValueVector(const std::vector<float>& vec)
{
	_isolinesVisualization.setTargetValueVector(vec);
}

inline float MainInterface::getIsolinesWidth()
{
	return _isolinesVisualization.getLineWidth();
}

inline void MainInterface::setIsolinesWidth(float width)
{
	_isolinesVisualization.setLineWidth(width);
}

inline float MainInterface::getHeightMax()
{
	return _heightPlotVisualization.getMaxHeight();
}

inline void MainInterface::setHeightMax(float height)
{
	_heightPlotVisualization.setMaxHeight(height);
}

inline float MainInterface::getGridHeightMax()
{
	return _heightPlotVisualization.getGridHeight();
}

inline void MainInterface::setGridHeightMax(float height)
{
	_heightPlotVisualization.setGridHeight(height);
}

inline bool MainInterface::getDrawHeightGrid()
{
	return _heightPlotVisualization.getDrawGrid();
}

inline void MainInterface::setDrawHeightGrid(bool value)
{
	_heightPlotVisualization.setDrawGrid(value);
}

inline bool MainInterface::getBlend()
{
	return _blend;
}

inline void MainInterface::setBlend(bool blend)
{
	_blend = blend;
}

inline void MainInterface::setHeightGridAlpha(float alpha)
{
	_heightPlotVisualization.setGridAlpha(alpha);
}

inline void MainInterface::setGranulomaThreshold(float refValue)
{
	_simulation->setAreaThreshold(refValue);
	_granulomaVisualization.setTargetValueVector(std::vector<float>(1, refValue));
}

inline int MainInterface::getTime()
{
	return _time;
}

inline int MainInterface::getSimTime()
{
	return _simTime;
}

inline void MainInterface::setTime(int time)
{
	_time = time;
}

inline void MainInterface::incTime(int dTime)
{
	_time += dTime;
}

inline void MainInterface::setSimTime(int simTime)
{
	_simTime = simTime;
}

#endif /* MAININTERFACE_H_ */
