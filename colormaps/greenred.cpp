/*
 * greenred.cpp
 *
 *  Created on: 17-dec-2008
 *      Author: S020751
 */

#include "greenred.h"
#include <float.h>

ColorMapGreenRed::ColorMapGreenRed()
{
}

ColorMapGreenRed::~ColorMapGreenRed()
{
}

void ColorMapGreenRed::map(float value, float& R, float& G, float& B) const
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
		h = 0.3333f;
		v = (0.5f - value) / 0.55f;
		s = 1.0f;
	}
	else
	{
		/* 0.5 - 1.0 red */
		h = 0.0f;
		v = (value - 0.5f) / 0.55f;
		s = 1.0f;
	}

	hsv2rgb(h, s, v, R, G, B);

	applyDeltas(R, G, B);
}
