/*
 * scalarintmtbdataset.h
 *
 *  Created on: Nov 24, 2009
 *      Author: mohammed
 */

#ifndef SCALARINTMTBDATASET_H_
#define SCALARINTMTBDATASET_H_

#include "scalardataset.h"

class ScalarIntMtbDataset : public ScalarDataset
{
public:
  ScalarIntMtbDataset();
  virtual ~ScalarIntMtbDataset();
  virtual float getScalar(const Simulation* pSimulation, int row, int col) const;
};

inline ScalarIntMtbDataset::ScalarIntMtbDataset()
  : ScalarDataset()
{
}

inline ScalarIntMtbDataset::~ScalarIntMtbDataset()
{
}

inline float ScalarIntMtbDataset::getScalar(const Simulation* pSimulation, int row, int col) const
{
  float sum = 0;
  for(size_t i=0; i<GrGrid::MAX_AGENTS_PER_CELL; i++)
    {
      const Mac* pMac = dynamic_cast<const Mac*>(pSimulation->getGrGrid().agent(row, col, i));
      sum += pMac ? pMac->getIntMtb() : 0;
    }
  return sum;
}

#endif /* SCALARINTMTBDATASET_H_ */
