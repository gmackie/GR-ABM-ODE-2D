/*
 * scalargrid.cpp
 *
 *  Created on: 21-sep-2008
 *      Author: S030858
 */

#include "scalargrid.h"
#include "simulation.h"
#include <fstream>
#include "scalardataset.h"

ScalarGrid::ScalarGrid(size_t _dim)
	: Grid(_dim)
	, _min(FLT_MAX)
	, _max(FLT_MIN)
	, _grid(_dim * _dim)
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

void ScalarGrid::evaluate(const Simulation* pSimulation, ScalarDataset* pScalarDataset, bool useNN)
{
	assert(pScalarDataset);

	_min = FLT_MAX;
	_max = FLT_MIN;

	for (size_t i = 0; i < _grid.size(); i++)
	{
		ScalarGridItem& item = _grid[i];
		if (useNN)
		{
			item.scalar = pScalarDataset->getScalarNN(pSimulation, item.pos);
		}
		else
		{
			item.scalar = pScalarDataset->getScalarBL(pSimulation, item.pos);
		}

		if (item.scalar > _max)
			_max = item.scalar;
		else if (item.scalar < _min)
			_min = item.scalar;
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
