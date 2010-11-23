/*
 * GrSimulationGrid.cpp
 *
 *  Created on: Nov 23, 2010
 *      Author: pwolberg
 */

#include "grsimulationgrid.h"

GrSimulationGrid::GrSimulationGrid()
	: _grid1()
	, _grid2()
	, _pCurrentGrid(&_grid1)
	, _pNextGrid(&_grid2)
{

}

GrSimulationGrid::~GrSimulationGrid() {

}
