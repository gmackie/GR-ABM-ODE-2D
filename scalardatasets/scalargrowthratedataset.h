/*
 * ScalarGrowthRateDataset.h
 *
 * A scalar dataset for determining a granuloma boundary by cell density.
 *
 * This uses an idea by Mohammad Fallahi
 *
 *  Created on: Apr 29, 2010
 *      Author: Paul Wolberg
 */

#ifndef SCALARGROWTHRATEDATASET_H_
#define SCALARGROWTHRATEDATASET_H_

#include "scalardataset.h"

class ScalarGrowthRateDataset: public ScalarDataset
{
public:
  ScalarGrowthRateDataset();
  virtual ~ScalarGrowthRateDataset();
  virtual float getScalar(const Simulation* pSimulation, int row, int col) const;
};

inline ScalarGrowthRateDataset::ScalarGrowthRateDataset()
  : ScalarDataset()
{
}

inline ScalarGrowthRateDataset::~ScalarGrowthRateDataset()
{
}

inline float ScalarGrowthRateDataset::getScalar(const Simulation* pSimulation, int row, int col) const
{

  // Check the micro-compartment's Moore neighborhood (including the micro-compartment itself).
  const GrGrid& g = pSimulation->getGrGrid();
  float sum=0.0;
  for (size_t k=0;k<GrGrid::MAX_AGENTS_PER_CELL;k++)
    if(Mac::isMac(g.agent(row, col, k)))
      sum += static_cast<const Mac*>(g.agent(row, col, k))->getGrowthRate();

  return sum;
}


#endif /* SCALARGROWTHRATEDATASET_H_ */
