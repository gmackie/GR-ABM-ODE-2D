/*
 * grgrid.h
 *
 *  Created on: 02-nov-2009
 *      Author: M. El-Kebir
 */

#ifndef GRGRID_H
#define GRGRID_H

#include "gr.h"
#include "gridcell.h"
#include "params.h"

class GrGrid
{
private:
	GridCell _grid[NROWS][NCOLS];
	GridCellPtrList _sources;

public:
	GrGrid();
	~GrGrid();
	GridCell& operator ()(int row, int col);
	GridCell operator ()(int row, int col) const;
	GridCellPtrList& getSources();
	void initSources();
	void serialize(std::ostream& out) const;
	void deserialize(std::istream& in);
};

inline GridCell& GrGrid::operator ()(int row, int col)
{
	assert(0 <= row && row < NROWS);
	assert(0 <= col && col < NCOLS);
	return _grid[row][col];
}

inline GridCell GrGrid::operator ()(int row, int col) const
{
	assert(0 <= row && row < NROWS);
	assert(0 <= col && col < NCOLS);
	return _grid[row][col];
}

inline GridCellPtrList& GrGrid::getSources()
{
	return _sources;
}

#endif /* GRID_H */
