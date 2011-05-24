/*
 * grid.cpp
 *
 *  Created on: 02-nov-2009
 *      Author: M. El-Kebir
 */

#include "grgrid.h"
#include <algorithm>

GrGrid::GrGrid()
	: _grid()
	, _sources()
{
	// reinitialize _grid
	for (int i = 0; i < NROWS; i++)
		for (int j = 0; j < NCOLS; j++)
			_grid[i][j].setRowCol(i, j);
}

GrGrid::~GrGrid()
{
}

void GrGrid::initSources()
{
	int nSources = _PARAM(PARAM_GR_NR_SOURCES);
	int n = (int)floor(sqrtf((float) nSources));

	int dRow = NROWS / n;
	int dCol = NCOLS / n;

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			int row = g_Rand.getInt(dRow);
			int col = g_Rand.getInt(dCol);

			GridCell* pGridCell = &_grid[dRow * i + row][dCol * j + col];

			assert(!pGridCell->isSource());

			pGridCell->setSource();
			_sources.push_back(pGridCell);
		}
	}

	// pick remaining sources
	while ((int)_sources.size() < nSources)
	{
		int row = g_Rand.getInt(NROWS);
		int col = g_Rand.getInt(NCOLS);

		GridCell* pGridCell = &_grid[row][col];
		if (!pGridCell->isSource())
		{
			pGridCell->setSource();
			_sources.push_back(pGridCell);
		}
	}
}

void GrGrid::shuffleSources()
{
	random_shuffle(_sources.begin(), _sources.end(), g_Rand);
}

void GrGrid::serialize(std::ostream& out) const
{
	assert(out.good());

	out << NROWS << std::endl;
	out << NCOLS << std::endl;

	for (int i = 0; i < NROWS; i++)
	{
		for (int j = 0; j < NCOLS; j++)
		{
			_grid[i][j].serialize(out);
		}
	}
}

void GrGrid::deserialize(std::istream& in)
{
	assert(in.good());

	int savedNrows, savedNcols;

	in >> savedNrows;
	in >> savedNcols;

	if (savedNrows != NROWS || savedNcols != NCOLS)
	{
		std::cerr << "Error deserializing a GrGrid:"<< std::endl;
		std::cerr << "The number of rows and cols in the saved grid, (" << savedNrows << ","  << savedNcols
				  << "), does not match the current grid (" << NROWS << ","  << NCOLS << ")."<< std::endl;
		exit(1);
	}

	_sources.clear();
	for (int i = 0; i < NROWS; i++)
	{
		for (int j = 0; j < NCOLS; j++)
		{
			_grid[i][j].deserialize(in);
			
			if (_grid[i][j].isSource())
				_sources.push_back(&_grid[i][j]);
		}
	}
}
