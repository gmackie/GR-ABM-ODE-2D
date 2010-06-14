/*
 * scalarextmtbdataset.h
 *
 *  Created on: Nov 24, 2009
 *      Author: mohammed
 */

#ifndef SCALAREXTMTBDATASET_H_
#define SCALAREXTMTBDATASET_H_

#include "scalardataset.h"

class ScalarExtMtbDataset : public ScalarDataset
{
public:
	ScalarExtMtbDataset();
	virtual ~ScalarExtMtbDataset();
	virtual float getScalar(const Simulation* pSimulation, int row, int col) const;
};

inline ScalarExtMtbDataset::ScalarExtMtbDataset()
	: ScalarDataset()
{
}

inline ScalarExtMtbDataset::~ScalarExtMtbDataset()
{
}

inline float ScalarExtMtbDataset::getScalar(const Simulation* pSimulation, int row, int col) const
{
	return pSimulation->getGrGrid()(row, col).getExtMtb();
}

#endif /* SCALAREXTMTBDATASET_H_ */
