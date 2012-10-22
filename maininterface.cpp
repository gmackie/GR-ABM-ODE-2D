/*
 * maininterface.cpp
 *
 *  Created on: 21-sep-2008
 *      Author: S030858
 */

#include "maininterface.h"
#include "scalardatasets/scalartnfdataset.h"
#include "scalardatasets/scalartnfattrextmtb.h"
#include "colormaps/fire.h"
#include "colormaps/blackwhite.h"
#include "colormaps/rainbow.h"
#include "colormaps/coolwarm.h"
#include "colormaps/fixed.h"
#include "scalardatasets/scalaragentgridbase.h"
#include <stdio.h>

MainInterface::MainInterface(const Pos& dim, Visualization* pAgentVisualization, ScalarAgentGridBase* pScalarAgentGrid)
  : _simulation(new Simulation(dim))
  , _scalarSmokeNormalizer(_MIN_CLAMP_VALUE, _MAX_CLAMP_VALUE)
  , _scalarGlyphNormalizer(_MIN_CLAMP_VALUE, _MAX_CLAMP_VALUE)
  , _scalarIsolinesNormalizer(_MIN_CLAMP_VALUE, _MAX_CLAMP_VALUE)
  , _scalarHeightNormalizer(_MIN_CLAMP_VALUE, _MAX_CLAMP_VALUE)
  , _scalarHeightColorNormalizer(_MIN_CLAMP_VALUE, _MAX_CLAMP_VALUE)
  , _scalarGranulomaNormalizer(_MIN_CLAMP_VALUE, _MAX_CLAMP_VALUE)
  , _pScalarSmokeDataset(new ScalarTnfDataset())
  , _pScalarGlyphDataset(new ScalarTnfDataset())
  , _pVectorGlyphDataset(NULL)
  , _pScalarHeightDataset(new ScalarTnfDataset())
  , _pScalarHeightColorDataset(new ScalarTnfDataset())
  , _pScalarGranulomaDataset(new ScalarTnfAttrExtMtb())
  , _pScalarAgentGrid(pScalarAgentGrid)
  , _scalarSmokeGrid(dim.x)
  , _scalarGlyphGrid(dim.x)
  , _vectorGlyphGrid(dim.x)
  , _scalarHeightGrid(dim.x)
  , _scalarHeightColorGrid(dim.x)
  , _scalarGranulomaGrid(dim.x)
  , _pAgentsVisualization(pAgentVisualization)
  , _invisibleQuadVisualization(dim.x)
  , _smokeVisualization(dim.x, &_scalarSmokeNormalizer, &_scalarSmokeGrid)
  , _glyphVisualization(dim.x, &_scalarGlyphNormalizer, &_scalarGlyphGrid, &_vectorGlyphGrid)
  , _isolinesVisualization(dim.x, &_scalarIsolinesNormalizer, &_scalarSmokeGrid)
  , _heightPlotVisualization(dim.x, &_scalarHeightNormalizer, &_scalarHeightColorNormalizer, &_scalarHeightGrid, &_scalarHeightColorGrid,
                             new VectorGradientDataset(_pScalarHeightDataset, true))
  , _granulomaVisualization(dim.x, &_scalarGranulomaNormalizer, &_scalarGranulomaGrid)
  , _pColorMapSmoke(new ColorMapRainbow())
  , _pColorMapGlyphs(new ColorMapBlackWhite())
  , _pColorMapIsolines(new ColorMapRainbow())
  , _pColorMapHeightPlot(new ColorMapRainbow())
  , _pColorMapGranuloma(new ColorMapFixed(_GR_BORDER_COLOR))
  , _drawAgents(_DRAW_AGENTS)
  , _drawSmoke(_DRAW_SMOKE)
  , _drawGlyphs(_DRAW_GLYPHS)
  , _drawIsolines(_DRAW_ISOLINES)
  , _drawHeightPlot(_DRAW_HEIGHTPLOT)
  , _drawGranuloma(_DRAW_GRANULOMA_BORDER)
  , _glyphGridNN(true)
  , _smokeGridNN(true)
  , _heightGridNN(true)
  , _blend(false)
  , _time(0)
  , _simTime(0)
{
  _pVectorGlyphDataset = new VectorGradientDataset(new ScalarTnfDataset(), getVectorGlyphGrid());

  /* setup granulomaVisualization */
  _granulomaVisualization.setTargetValueVector(std::vector<float>(1, _AREA_THRESHOLD));
  _granulomaVisualization.setLineWidth(_GR_BORDER_SIZE);
  _pColorMapGranuloma->setAlpha(_GR_BORDER_ALPHA);
}

MainInterface::~MainInterface()
{
  delete _simulation;
  delete _pScalarSmokeDataset;
  delete _pScalarGlyphDataset;
  delete _pVectorGlyphDataset;
  delete _pScalarHeightColorDataset;
  delete _pScalarHeightDataset;
  delete _pScalarGranulomaDataset;
  delete _pColorMapSmoke;
  delete _pColorMapGlyphs;
  delete _pColorMapIsolines;
  delete _pColorMapHeightPlot;
  delete _pColorMapGranuloma;
}

void MainInterface::visualize()
{
  if (_drawGlyphs)
    {
      _glyphVisualization.visualize(_blend, NULL, _pColorMapGlyphs);
    }

  if (_drawIsolines)
    {
      _isolinesVisualization.visualize(_blend, NULL, _pColorMapIsolines);
    }

  if (_drawGranuloma)
    {
      _granulomaVisualization.visualize(_blend, NULL, _pColorMapGranuloma);
    }

  if (_drawAgents)
    {
      _pAgentsVisualization->visualize(_blend, NULL, NULL);
    }

  if (_drawSmoke)
    {
      _smokeVisualization.visualize(_blend, NULL, _pColorMapSmoke);
    }
  else if (_drawHeightPlot)
    {
      _heightPlotVisualization.visualize(_blend, _simulation, _pColorMapHeightPlot);
    }
  else
    {
      _invisibleQuadVisualization.visualize(_blend, NULL, NULL);
    }
}

void MainInterface::updateGrids()
{

  if (_drawAgents)
    {
      _pScalarAgentGrid->evaluate(_simulation);
    }

  _scalarGranulomaGrid.evaluate(_simulation, _pScalarGranulomaDataset, true);

  if (_drawSmoke || _drawIsolines)
    {
      _scalarSmokeGrid.evaluate(_simulation, _pScalarSmokeDataset, _smokeGridNN);
    }

  if (_drawGlyphs)
    {
      _scalarGlyphGrid.evaluate(_simulation, _pScalarGlyphDataset, _glyphGridNN);
      _vectorGlyphGrid.evaluate(_simulation, _pVectorGlyphDataset, _glyphGridNN);
    }

  if (_drawHeightPlot)
    {
      _scalarHeightGrid.evaluate(_simulation, _pScalarHeightDataset, _heightGridNN);
      _scalarHeightColorGrid.evaluate(_simulation, _pScalarHeightColorDataset, _heightGridNN);
    }

  if (!_scalarSmokeNormalizer.getClamping())
    {
      _scalarSmokeNormalizer.setMin(_scalarSmokeGrid.getMin());
      _scalarSmokeNormalizer.setMax(_scalarSmokeGrid.getMax());
    }

  if (!_scalarGlyphNormalizer.getClamping())
    {
      _scalarGlyphNormalizer.setMin(_scalarGlyphGrid.getMin());
      _scalarGlyphNormalizer.setMax(_scalarGlyphGrid.getMax());
    }

  if (!_scalarIsolinesNormalizer.getClamping())
    {
      _scalarIsolinesNormalizer.setMin(_scalarSmokeGrid.getMin());
      _scalarIsolinesNormalizer.setMax(_scalarSmokeGrid.getMax());
    }

  if (!_scalarHeightNormalizer.getClamping())
    {
      _scalarHeightNormalizer.setMin(_scalarHeightGrid.getMin());
      _scalarHeightNormalizer.setMax(_scalarHeightGrid.getMax());
    }

  if (!_scalarHeightColorNormalizer.getClamping())
    {
      _scalarHeightColorNormalizer.setMin(_scalarHeightColorGrid.getMin());
      _scalarHeightColorNormalizer.setMax(_scalarHeightColorGrid.getMax());
    }

}

void MainInterface::serializeGrids(std::ofstream& outFile)
{
  outFile << "Scalar Smoke Grid" << std::endl;
  _scalarSmokeGrid.serialize(outFile);

  outFile << "Scalar Glyph Grid" << std::endl;
  _scalarGlyphGrid.serialize(outFile);

  outFile << "Vector Glyph Grid" << std::endl;
  _vectorGlyphGrid.serialize(outFile);

  outFile << "Scalar Height Grid" << std::endl;
  _scalarHeightGrid.serialize(outFile);

  outFile << "Scalar Height Color Grid" << std::endl;
  _scalarHeightColorGrid.serialize(outFile);
}
