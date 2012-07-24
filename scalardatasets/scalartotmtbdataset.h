/*
 * scalarextmtbdataset.h
 *
 *  Created on: Nov 24, 2009
 *      Author: mohammed
 */

#ifndef SCALARTOTMTBDATASET_H_
#define SCALARTOTMTBDATASET_H_

#include "scalardataset.h"

class ScalarTotMtbDataset : public ScalarDataset
{
public:
  ScalarTotMtbDataset();
  virtual ~ScalarTotMtbDataset();
  virtual float getScalar(const Simulation* pSimulation, int row, int col) const;
};

inline ScalarTotMtbDataset::ScalarTotMtbDataset()
  : ScalarDataset()
{
}

inline ScalarTotMtbDataset::~ScalarTotMtbDataset()
{
}

inline float ScalarTotMtbDataset::getScalar(const Simulation* pSimulation, int row, int col) const
{
  float totMtb = 0;
  for(size_t i=0; i<GrGrid::MAX_AGENTS_PER_CELL; i++)
    {
      const Mac* pMac = dynamic_cast<const Mac*>(pSimulation->getGrGrid().agent(row, col, i));
      totMtb += pMac ? pMac->getIntMtb() : 0;
    }

  totMtb += pSimulation->getGrGrid().extMTB(row, col);

  return totMtb;
}

#endif /* SCALAREXTMTBDATASET_H_ */
