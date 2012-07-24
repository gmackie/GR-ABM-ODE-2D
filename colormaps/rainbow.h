/*
 * rainbow.h
 *
 *  Created on: 8-sep-2008
 *      Author: s030858
 */

#ifndef RAINBOW_H_
#define RAINBOW_H_

#include "colormap.h"
#include <QString>

class ColorMapRainbow : public ColorMap
{
public:
  ColorMapRainbow();
  virtual ~ColorMapRainbow();
  virtual void map(float value, float& R, float& G, float& B) const;
  virtual QString getName() const;
};

inline QString ColorMapRainbow::getName() const
{
  return QString("Rainbow");
}

#endif /* RAINBOW_H_ */
