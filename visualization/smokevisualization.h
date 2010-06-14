/*
 * smokevisualization.h
 *
 *  Created on: 21-sep-2008
 *      Author: S030858
 */

#ifndef SMOKEVISUALIZATION_H_
#define SMOKEVISUALIZATION_H_

#include "colormaps/colormap.h"
#include "scalardatasets/scalargrid.h"
#include "scalardatasets/scalarnormalizer.h"
#include "visualization.h"

class SmokeVisualization : public Visualization
{
public:
	SmokeVisualization(int DIM, const ScalarNormalizer* pScalarNormalizer, const ScalarGrid* pScalarGrid);
	virtual ~SmokeVisualization();
	void visualize(bool blend, const Simulation* pSimulation, const ColorMap* pColorMap) const;

private:
	const ScalarNormalizer* _pScalarNormalizer;
	const ScalarGrid* _pScalarGrid;
};

#endif /* SMOKEVISUALIZATION_H_ */
