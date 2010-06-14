/*
 * scalarattractantdataset.h
 *
 *  Created on: Apr 19, 2010
 *      Author: mohammed
 */

#ifndef SCALARATTRACTANTDATASET_H_
#define SCALARATTRACTANTDATASET_H_

class ScalarAttractantDataset : public ScalarDataset
{
public:
	ScalarAttractantDataset();
	virtual ~ScalarAttractantDataset();
	virtual float getScalar(const Simulation* pSimulation, int row, int col) const;
};

inline ScalarAttractantDataset::ScalarAttractantDataset()
	: ScalarDataset()
{
}

inline ScalarAttractantDataset::~ScalarAttractantDataset()
{
}

inline float ScalarAttractantDataset::getScalar(const Simulation* pSimulation, int row, int col) const
{
	return pSimulation->getGrGrid()(row, col).getMacAttractant();
}


#endif /* SCALARATTRACTANTDATASET_H_ */
