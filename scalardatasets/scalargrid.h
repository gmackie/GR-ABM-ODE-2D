/*
 * scalargrid.h
 *
 *  Created on: 21-sep-2008
 *      Author: S030858
 */

#ifndef SCALARGRID_H_
#define SCALARGRID_H_

#include "grviz.h"
#include "simulation.h"
#include "datasets/grid.h"
#include "vectordatasets/vector.h"
#include <fstream>

class ScalarDataset;

struct ScalarGridItem
{
  vec2f pos;
  float scalar;
};

class ScalarGrid : public Grid
{
private:
  float _min;
  float _max;
  std::vector<ScalarGridItem> _grid;

public:
  ScalarGrid(size_t dim);
  ~ScalarGrid();
  void evaluate(const Simulation* pSimulation, ScalarDataset* pScalarDataset, bool useNN);
  int count() const;
  float getMin() const;
  float getMax() const;
  const std::vector<ScalarGridItem>& getGrid() const;
  void serialize(std::ofstream& outFile) const;
};

inline int ScalarGrid::count() const
{
  return (int)_grid.size();
}

inline float ScalarGrid::getMin() const
{
  return _min;
}

inline float ScalarGrid::getMax() const
{
  return _max;
}

inline const std::vector<ScalarGridItem>& ScalarGrid::getGrid() const
{
  return _grid;
}

#endif /* SCALARGRID_H_ */
