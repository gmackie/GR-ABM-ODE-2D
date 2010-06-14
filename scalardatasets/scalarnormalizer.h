/*
 * scalarnormalizer.h
 *
 *  Created on: 12-sep-2008
 *      Author: s030858
 */

#ifndef SCALARNORMALIZER_H_
#define SCALARNORMALIZER_H_

#include "scalardataset.h"

class ScalarNormalizer
{
private:
	float _min;
	float _max;
	bool _clamping;

public:
	ScalarNormalizer(float min, float max);
	float normalize(float value) const;
	float denormalize(float normalizedValue) const;
	float getMin() const;
	float getMax() const;
	bool getClamping() const;
	void setMin(float f);
	void setMax(float f);
	void setClamping(bool clamping);
	virtual ~ScalarNormalizer();
};

inline float ScalarNormalizer::getMin() const
{
	return _min;
}

inline float ScalarNormalizer::getMax() const
{
	return _max;
}

inline bool ScalarNormalizer::getClamping() const
{
	return _clamping;
}

inline void ScalarNormalizer::setMin(float f)
{
	_min = f;
}

inline void ScalarNormalizer::setMax(float f)
{
	_max = f;
}

inline void ScalarNormalizer::setClamping(bool clamping)
{
	_clamping = clamping;
}

#endif /* SCALARNORMALIZER_H_ */
