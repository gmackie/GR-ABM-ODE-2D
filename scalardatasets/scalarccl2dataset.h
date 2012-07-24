/*
 * scalarccl2dataset.h
 *
 *  Created on: Nov 24, 2009
 *      Author: mohammed
 */

#ifndef SCALARCCL2DATASET_H_
#define SCALARCCL2DATASET_H_

#include "scalardataset.h"

class ScalarCcl2Dataset : public ScalarDataset
{
public:
  ScalarCcl2Dataset();
  virtual ~ScalarCcl2Dataset();
  virtual float getScalar(const Simulation* pSimulation, int row, int col) const;
};

inline ScalarCcl2Dataset::ScalarCcl2Dataset()
  : ScalarDataset()
{
}

inline ScalarCcl2Dataset::~ScalarCcl2Dataset()
{
}

inline float ScalarCcl2Dataset::getScalar(const Simulation* pSimulation, int row, int col) const
{
  return pSimulation->getGrGrid().CCL2(row, col);
}

#endif /* SCALARCCL2DATASET_H_ */
