/*
 * scalaril10dataset.h
 *
 *  Created on: Mar 7, 2012
 *      Author: mohammed
 */

#ifndef SCALARIL10DATASET_H_
#define SCALARIL10DATASET_H_

#include "scalardataset.h"

class ScalarIl10Dataset : public ScalarDataset
{
public:
	ScalarIl10Dataset();
	virtual ~ScalarIl10Dataset();
	virtual float getScalar(const Simulation* pSimulation, int row, int col) const;
};

inline ScalarIl10Dataset::ScalarIl10Dataset()
	: ScalarDataset()
{
}

inline ScalarIl10Dataset::~ScalarIl10Dataset()
{
}

inline float ScalarIl10Dataset::getScalar(const Simulation* pSimulation, int row, int col) const
{
	return pSimulation->getGrGrid().il10(row, col);
}

#endif /* scalarIl10DATASET_H_ */
