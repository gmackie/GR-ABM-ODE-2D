/*
 * blackwhite.cpp
 *
 *  Created on: 8-sep-2008
 *      Author: s030858
 */

#include "blackwhite.h"
#include <assert.h>
#include <float.h>

ColorMapBlackWhite::ColorMapBlackWhite()
{
}

ColorMapBlackWhite::~ColorMapBlackWhite()
{
}

void ColorMapBlackWhite::map(float value, float& R, float& G, float& B) const
{
  if (_invert)
    {
      invert(value);
    }
  bin(value);

  R = G = B = value;

  applyDeltas(R, G, B);
}
