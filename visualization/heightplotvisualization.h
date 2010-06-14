/*
 * heightplotvisualization.h
 *
 *  Created on: 18-nov-2008
 *      Author: s030858
 */

#ifndef HEIGHTPLOTVISUALIZATION_H_
#define HEIGHTPLOTVISUALIZATION_H_

#include "colormaps/colormap.h"
#include "scalardatasets/scalargrid.h"
#include "scalardatasets/scalarnormalizer.h"
#include "vectordatasets/vectorgradientdataset.h"
#include "visualization.h"

class HeightPlotVisualization : public Visualization
{
public:
	HeightPlotVisualization(int DIM, const ScalarNormalizer* pScalarHeightNormalizer,
			const ScalarNormalizer* pScalarHeightColorNormalizer,
			const ScalarGrid* pScalarHeightGrid, const ScalarGrid* pScalarHeightColorGrid,
			VectorGradientDataset* pGradientDataset);
	virtual ~HeightPlotVisualization();
	void visualize(bool blend, const Simulation* pSimulation, const ColorMap* pColorMap) const;
	float getMaxHeight() const;
	void setMaxHeight(float value);
	float getGridHeight() const;
	void setGridHeight(float value);
	bool getDrawGrid() const;
	void setDrawGrid(bool value);
	void setGradientDataset(VectorGradientDataset* pGradientDataset);
	float getGridAlpha() const;
	void setGridAlpha(float alpha);

private:
	const ScalarNormalizer* _pScalarHeightNormalizer;
	const ScalarNormalizer* _pScalarHeightColorNormalizer;
	const ScalarGrid* _pScalarHeightGrid;
	const ScalarGrid* _pScalarHeightColorGrid;
	float _maxHeight;
	bool _drawGrid;
	float _gridHeight;
	float _gridAlpha;
	VectorGradientDataset* _pGradientDataset;
	void drawGrid() const;
};

inline float HeightPlotVisualization::getGridAlpha() const
{
	return _gridAlpha;
}

inline void HeightPlotVisualization::setGridAlpha(float alpha)
{
	_gridAlpha = alpha;
}

inline float HeightPlotVisualization::getMaxHeight() const
{
	return _maxHeight;
}

inline void HeightPlotVisualization::setMaxHeight(float value)
{
	_maxHeight = value;
}

inline float HeightPlotVisualization::getGridHeight() const
{
	return _gridHeight;
}

inline void HeightPlotVisualization::setGridHeight(float value)
{
	_gridHeight = value;
}

inline bool HeightPlotVisualization::getDrawGrid() const
{
	return _drawGrid;
}

inline void HeightPlotVisualization::setDrawGrid(bool value)
{
	_drawGrid = value;
}

inline void HeightPlotVisualization::setGradientDataset(VectorGradientDataset* pGradientDataset)
{
	delete _pGradientDataset;
	_pGradientDataset = pGradientDataset;
}

#endif /* HEIGHTPLOTVISUALIZATION_H_ */
