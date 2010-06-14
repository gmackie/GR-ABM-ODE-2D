/*
 * fire.cpp
 *
 *  Created on: 12-sep-2008
 *      Author: s030858
 */

#include "fire.h"
#include <algorithm>
#include <float.h>

ColorMapFire::ColorMapFire()
{
}

ColorMapFire::~ColorMapFire()
{
}

void ColorMapFire::map(float value, float& R, float& G, float& B) const
{
	if (_invert)
	{
		invert(value);
	}
	bin(value);

	float magic = 85.0f / 255.0f;
	R = std::min(value, magic) * 3.0f;
	value = std::max(value - magic, 0.0f);
	G = std::min(value, magic) * 3.0f;
	value = std::max(value - magic, 0.0f);
	B = std::min(value, magic) * 3.0f;

	applyDeltas(R, G, B);
}
