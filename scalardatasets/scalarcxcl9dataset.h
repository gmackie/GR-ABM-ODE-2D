/*
 * scalarcxcl9dataset.h
 *
 *  Created on: Nov 24, 2009
 *      Author: mohammed
 */

#ifndef SCALARCXCL9DATASET_H_
#define SCALARCXCL9DATASET_H_

#include "scalardataset.h"

class ScalarCxcl9Dataset : public ScalarDataset
{
public:
	ScalarCxcl9Dataset();
	virtual ~ScalarCxcl9Dataset();
	virtual float getScalar(const Simulation* pSimulation, int row, int col) const;
};

inline ScalarCxcl9Dataset::ScalarCxcl9Dataset()
	: ScalarDataset()
{
}

inline ScalarCxcl9Dataset::~ScalarCxcl9Dataset()
{
}

inline float ScalarCxcl9Dataset::getScalar(const Simulation* pSimulation, int row, int col) const
{
	return pSimulation->getGrGrid().CXCL9(row, col);
}

#endif /* SCALARCXCL9DATASET_H_ */
