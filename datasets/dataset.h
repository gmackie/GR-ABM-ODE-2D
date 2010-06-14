/*
 * dataset.h
 *
 *  Created on: 2-nov-2008
 *      Author: s030858
 */

#ifndef DATASET_H_
#define DATASET_H_

#include <fstream>

class Dataset
{
protected:
	int moduloDIM(int value) const;
	float moduloDIM(float value) const;

public:
	Dataset();
	virtual ~Dataset();

};

inline Dataset::Dataset()
{
}

inline Dataset::~Dataset()
{
}

inline int Dataset::moduloDIM(int value) const
{
	if (value < 0)
	{
		value += (value * -1 / Simulation::_DIM + 1) * Simulation::_DIM;
		return value;
	}
	else if (value >= Simulation::_DIM)
	{
		return value % Simulation::_DIM;
	}
	else
	{
		return value;
	}
}

inline float Dataset::moduloDIM(float value) const
{
	while (value < 0.0f)
	{
		value += Simulation::_DIM;
	}
	while (value >= Simulation::_DIM)
	{
		value -= Simulation::_DIM;
	}

	return value;
}

#endif /* DATASET_H_ */
