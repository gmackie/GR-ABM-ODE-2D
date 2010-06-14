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
	const GridCell& cell = pSimulation->getGrGrid()(row, col);

	const Mac* pMac0 = dynamic_cast<const Mac*>(cell.getAgent(0));
	const Mac* pMac1 = dynamic_cast<const Mac*>(cell.getAgent(1));
	return (pMac0 ? pMac0->getIntMtb() : 0) + (pMac1 ? pMac1->getIntMtb() : 0);
}

#endif /* SCALARINTMTBDATASET_H_ */
