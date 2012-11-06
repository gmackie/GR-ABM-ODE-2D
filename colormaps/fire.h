/*
 * fire.h
 *
 *  Created on: 12-sep-2008
 *      Author: s030858
 */

#ifndef FIRE_H_
#define FIRE_H_

#include "colormap.h"
#include <QString>

class ColorMapFire : public ColorMap
{
public:
  ColorMapFire();
  ~ColorMapFire();
  void map(float value, float& R, float& G, float& B) const;
  QString getName() const;
  COLORMAP getType() const;
};

inline QString ColorMapFire::getName() const
{
  return QString("Fire");
}

inline COLORMAP ColorMapFire::getType() const
{
  return CMAPFIRE;
}

#endif /* FIRE_H_ */
