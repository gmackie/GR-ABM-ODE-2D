/*
 * scalarnormalizer.cpp
 *
 *  Created on: 12-sep-2008
 *      Author: s030858
 */

#include "scalarnormalizer.h"
#include <assert.h>

ScalarNormalizer::ScalarNormalizer(float min, float max)
	: _min(min)
	, _max(max)
	, _clamping(true)
{
}

ScalarNormalizer::~ScalarNormalizer()
{
}

float ScalarNormalizer::normalize(float value) const
{
	if (value < _min)
		value = _min;
	else if (value > _max)
		value = _max;

	return (value - _min) / (_max - _min);
}

float ScalarNormalizer::denormalize(float normalizedValue) const
{
	assert(0 <= normalizedValue && normalizedValue <= 1);
	return normalizedValue * (_max - _min) + _min;
}
