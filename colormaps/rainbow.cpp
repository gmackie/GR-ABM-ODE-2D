/*
 * rainbow.cpp
 *
 *  Created on: 8-sep-2008
 *      Author: s030858
 */

#include "rainbow.h"

ColorMapRainbow::ColorMapRainbow()
{
}

ColorMapRainbow::~ColorMapRainbow()
{
}

void ColorMapRainbow::map(float value, float& R, float& G, float& B) const
{
	if (_invert)
	{
		invert(value);
	}
	bin(value);

	const float dx = 0.8f;
	if (value < 0)
		value = 0;
	if (value > 1)
		value = 1;
	value = (6 - 2 * dx) * value + dx;
	R = std::max(0.0f, (3 - fabsf(value - 4) - fabsf(value - 5)) / 2.0f);
	G = std::max(0.0f, (4 - fabsf(value - 2) - fabsf(value - 4)) / 2.0f);
	B = std::max(0.0f, (3 - fabsf(value - 1) - fabsf(value - 2)) / 2.0f);

	applyDeltas(R, G, B);
}
