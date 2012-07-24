/*
 * grid.h
 *
 *  Created on: 2-nov-2008
 *      Author: s030858
 */

#ifndef GRID_H_
#define GRID_H_

#include "simulation.h"

class Grid
{
protected:
  const int _DIM;

public:
  Grid(const int DIM);
  virtual ~Grid();
};

inline Grid::Grid(const int DIM)
  : _DIM(DIM)
{
}

inline Grid::~Grid()
{
}

#endif /* GRID_H_ */
