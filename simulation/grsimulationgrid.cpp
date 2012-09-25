/*
 * GrSimulationGrid.cpp
 *
 *  Created on: Nov 23, 2010
 *      Author: pwolberg
 */

#include "grsimulationgrid.h"

GrSimulationGrid::GrSimulationGrid(const Pos& dim)
{
  _pCurrentGrid = (new GrGrid(dim));
  _pNextGrid = (new GrGrid(dim));

}

GrSimulationGrid::~GrSimulationGrid()
{
  delete _pCurrentGrid;
  delete _pNextGrid;
}
