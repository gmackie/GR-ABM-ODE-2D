/*
 * GrSimulationGrid.cpp
 *
 *  Created on: Nov 23, 2010
 *      Author: pwolberg
 */

#include "grsimulationgrid.h"

GrSimulationGrid::GrSimulationGrid()
	: _pCurrentGrid(new GrGrid())
	, _pNextGrid(new GrGrid())
{

}

GrSimulationGrid::~GrSimulationGrid() {

}
