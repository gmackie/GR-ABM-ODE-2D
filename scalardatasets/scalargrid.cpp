/*
 * scalargrid.cpp
 *
 *  Created on: 21-sep-2008
 *      Author: S030858
 */

#include "scalargrid.h"
#include "simulation.h"
#include <fstream>

ScalarGrid::ScalarGrid()
	: Grid(Simulation::_DIM)
	, _min(FLT_MAX)
	, _max(FLT_MIN)
	, _grid(_DIM * _DIM)
{
	for (int i = 0; i < _DIM; i++) // row
	{
		for (int j = 0; j < _DIM; j++) // col
		{
			ScalarGridItem& item = _grid[j + i * _DIM];
			item.pos[0] = j;
			item.pos[1] = i;
			item.scalar = 0.0;
		}
	}
}

void ScalarGrid::serialize(std::ofstream& outFile) const
{
	outFile << _DIM << ";" << _DIM << std::endl;

	for (int i = 0; i < _DIM; i++) // row
	{
		for (int j = 0; j < _DIM; j++) // col
		{
			if (j != 0)
			{
				outFile << ";";
			}
			const ScalarGridItem& item = _grid[i * _DIM + j];
			outFile << item.scalar;
		}
		outFile << std::endl;
	}

	outFile << std::endl;
}

ScalarGrid::~ScalarGrid()
{
}
