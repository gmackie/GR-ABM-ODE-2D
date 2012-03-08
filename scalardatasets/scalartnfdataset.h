/*
 * scalartnfdataset.h
 *
 *  Created on: Nov 24, 2009
 *      Author: mohammed
 */

#ifndef SCALARTNFDATASET_H_
#define SCALARTNFDATASET_H_

#include "scalardataset.h"

class ScalarTnfDataset : public ScalarDataset
{
public:
	ScalarTnfDataset();
	virtual ~ScalarTnfDataset();
	virtual float getScalar(const Simulation* pSimulation, int row, int col) const;
};

inline ScalarTnfDataset::ScalarTnfDataset()
	: ScalarDataset()
{
}

inline ScalarTnfDataset::~ScalarTnfDataset()
{
}

inline float ScalarTnfDataset::getScalar(const Simulation* pSimulation, int row, int col) const
{
	return pSimulation->getGrGrid().TNF(row, col);
}

#endif /* SCALARTNFDATASET_H_ */
