/*
 * fixed.cpp
 *
 *  Created on: Jan 19, 2010
 *      Author: mohammed
 */

#include "fixed.h"

ColorMapFixed::ColorMapFixed(const QColor& color)
	: _color(color)
{
}

ColorMapFixed::~ColorMapFixed()
{
}

void ColorMapFixed::map(float, float& R, float& G, float& B) const
{
	qreal cpyR, cpyG, cpyB;
	_color.getRgbF(&cpyR, &cpyG, &cpyB);

	R = (float) cpyR;
	G = (float) cpyG;
	B = (float) cpyB;
}
