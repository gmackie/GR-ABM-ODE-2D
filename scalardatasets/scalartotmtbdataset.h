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
	const GridCell& cell = pSimulation->getGrGrid()(row, col);

	const Mac* pMac0 = dynamic_cast<const Mac*>(cell.getAgent(0));
	const Mac* pMac1 = dynamic_cast<const Mac*>(cell.getAgent(1));
	float intMtb =  (pMac0 ? pMac0->getIntMtb() : 0) + (pMac1 ? pMac1->getIntMtb() : 0);

	float totMtb = pSimulation->getGrGrid()(row, col).getExtMtb() + intMtb;

	return totMtb;
}

#endif /* SCALAREXTMTBDATASET_H_ */
