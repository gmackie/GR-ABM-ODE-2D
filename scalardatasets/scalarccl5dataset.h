/*
 * scalarccl5dataset.h
 *
 *  Created on: Nov 24, 2009
 *      Author: mohammed
 */

#ifndef SCALARCCL5DATASET_H_
#define SCALARCCL5DATASET_H_

#include "scalardataset.h"

class ScalarCcl5Dataset : public ScalarDataset
{
public:
	ScalarCcl5Dataset();
	virtual ~ScalarCcl5Dataset();
	virtual float getScalar(const Simulation* pSimulation, int row, int col) const;
};

inline ScalarCcl5Dataset::ScalarCcl5Dataset()
	: ScalarDataset()
{
}

inline ScalarCcl5Dataset::~ScalarCcl5Dataset()
{
}

inline float ScalarCcl5Dataset::getScalar(const Simulation* pSimulation, int row, int col) const
{
	return pSimulation->getGrGrid().CCL5(row, col);
}

#endif /* SCALARCCL5DATASET_H_ */
