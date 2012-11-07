/*
 * scalardataset.h
 *
 *  Created on: 12-sep-2008
 *      Author: s030858
 */

#ifndef SCALARDATASET_H_
#define SCALARDATASET_H_

#include "grviz.h"
#include "simulation.h"
#include "datasets/dataset.h"
#include "vectordatasets/vector.h"

class ScalarDataset : public Dataset
{
public:
  ScalarDataset();
  virtual ~ScalarDataset();
  virtual float getScalar(const Simulation* pSimulation, int row, int col) const = 0;
  virtual float getScalarNN(const Simulation* pSimulation, float x, float y) const;
  virtual float getScalarNN(const Simulation* pSimulation, const vec2f& pos) const;
  virtual float getScalarBL(const Simulation* pSimulation, float x, float y) const;
  virtual float getScalarBL(const Simulation* pSimulation, const vec2f& pos) const;
};

inline ScalarDataset::ScalarDataset()
  : Dataset()
{
}

inline ScalarDataset::~ScalarDataset()
{
}

inline float ScalarDataset::getScalarNN(const Simulation* pSimulation, float x, float y) const
{
  const Pos& dim = pSimulation->getSize();
  int col = moduloDIM((int)floor(x+0.5f), dim.x);
  int row = moduloDIM((int)floor(y+0.5f), dim.y);
  return getScalar(pSimulation, row, col);
}

inline float ScalarDataset::getScalarNN(const Simulation* pSimulation, const vec2f& pos) const
{
  return getScalarNN(pSimulation, pos[0], pos[1]);
}

inline float ScalarDataset::getScalarBL(const Simulation* pSimulation, float x, float y) const
{
  const Pos& dim = pSimulation->getSize();
  x = moduloDIM(x, dim.x);
  y = moduloDIM(y, dim.y);

  int col = (int)x;
  int row = (int)y;

  // make sure that we do not index out of bounds
  int sucRow = moduloDIM(row + 1, dim.x);
  int sucCol = moduloDIM(col + 1, dim.y);

  // first get the four scalar values of the cell
  float v1 = getScalar(pSimulation, row, col);
  float v2 = getScalar(pSimulation, row, sucCol);
  float v3 = getScalar(pSimulation, sucRow, col);
  float v4 = getScalar(pSimulation, sucRow, sucCol);

  // interpolate in the x-direction
  float t = x - (int)x;

  float q1;
  float q2;

  q1 = v1 + t*(v2 - v1);
  q2 = v3 + t*(v4 - v3);

  // interpolate in the y-direction
  t = y - (int)y;

  return q1 + t*(q2 - q1);
}

inline float ScalarDataset::getScalarBL(const Simulation* pSimulation, const vec2f& pos) const
{
  return getScalarBL(pSimulation, pos[0], pos[1]);
}

#endif /* SCALARDATASET_H_ */
