/*
 * coolwarm.cpp
 *
 *  Created on: 12-sep-2008
 *      Author: s030858
 */

#include "coolwarm.h"
#include <float.h>

ColorMapCoolWarm::ColorMapCoolWarm()
{
}

ColorMapCoolWarm::~ColorMapCoolWarm()
{
}

void ColorMapCoolWarm::map(float value, float& R, float& G, float& B) const
{
	if (_invert)
	{
		invert(value);
	}
	bin(value);

	float h, s, v;

	if (value < 0.5)
	{
		/* 0.0 - 0.5 blue*/
		h = 0.6667f;
		v = 1.0f;
		s = (0.5f - value) / 0.55f;
	}
	else
	{
		/* 0.5 - 1.0 red */
		h = 0.0f;
		v = 1.0f;
		s = (value - 0.5f) / 0.55f;
	}

	hsv2rgb(h, s, v, R, G, B);

	applyDeltas(R, G, B);
}
