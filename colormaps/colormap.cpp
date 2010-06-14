/*
 * colormap.cpp
 *
 *  Created on: 8-sep-2008
 *      Author: s030858
 */

#include "colormap.h"
#include <algorithm>
#include <float.h>

ColorMap::ColorMap()
	: _invert(false)
	, _nrBands(_COLORMAP_NR_BANDS)
	, _hueDelta(_COLORMAP_HUE_DELTA)
	, _satDelta(_COLORMAP_SAT_DELTA)
	, _valDelta(_COLORMAP_VAL_DELTA)
	, _alpha(_COLORMAP_ALPHA)

{
}

ColorMap::~ColorMap()
{
}

void ColorMap::bin(float& value) const
{
	value *= _nrBands;
	value = (int)value;

	if (value == _nrBands)
		value--;

	value /= _nrBands;
}

void ColorMap::applyDeltas(float& r, float& g, float& b) const
{
	float h, s, v;
	rgb2hsv(r, g, b, h, s, v);

	h = h + _hueDelta;
	while (h < 0)
	{
		h += 1;
	}
	while (h > 1)
	{
		h -= 1;
	}

	s = std::min(std::max(s + _satDelta, 0.0f), 1.0f);
	v = std::min(std::max(v + _valDelta, 0.0f), 1.0f);

	hsv2rgb(h, s, v, r, g, b);
}

void ColorMap::rgb2hsv(float r, float g, float b, float& h, float& s, float& v) const
{
	float M = std::max(r, std::max(g,b));
	float m = std::min(r, std::min(g,b));
	float d = M - m;

	v = M; 							// value = max(r, g, b)
	s = (M > 0.00001f) ? d / M : 0;	// saturation

	if (s == 0)
		h = 0.0f;					// achromatic case, hue = 0 by convention
	else
	{
		if (r == M)
			h = (g - b) / d;
		else if (g == M)
			h = 2 + (b - r) / d;
		else
			h = 4 + (r - g) / d;

		h /= 6;
		if (h < 0)
			h += 1;
	}
}

void ColorMap::hsv2rgb(float h, float s, float v, float& r, float& g, float& b) const
{
	int hueCase = (int)(h*6);
	float frac = 6 * h - hueCase;
	float lx = v * (1 - s);
	float ly = v * (1 - s * frac);
	float lz = v * (1 - s * (1 - frac));
	switch (hueCase)
	{
		case 0:
		case 6: // 0 < hue < 1/6
			r = v;
			g = lz;
			b = lx;
			break;
		case 1: // 1/6 < hue < 2/6
			r = ly;
			g = v;
			b = lx;
			break;
		case 2: // 2/6 < hue < 3/6
			r = lx;
			g = v;
			b = lz;
			break;
		case 3: // 3/6 < hue < 4/6
			r = lx;
			g = ly;
			b = v;
			break;
		case 4: // 4/6 < hue < 5/6
			r = lz;
			g = lx;
			b = v;
			break;
		case 5: // 5/6 < hue < 6/6
			r = v;
			g = lx;
			b = ly;
			break;
	}
}
