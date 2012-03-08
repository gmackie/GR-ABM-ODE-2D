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
	int moduloDIM(int value, int dim) const;
	float moduloDIM(float value, int dim) const;

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

inline int Dataset::moduloDIM(int value, int dim) const
{
	if (value < 0)
	{
		value += (value * -1 / dim + 1) * dim;
		return value;
	}
	else if (value >= dim)
	{
		return value % dim;
	}
	else
	{
		return value;
	}
}

inline float Dataset::moduloDIM(float value, int dim) const
{
	while (value < 0.0f)
	{
		value += dim;
	}
	while (value >= dim)
	{
		value -= dim;
	}

	return value;
}

#endif /* DATASET_H_ */
