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
  virtual ~ColorMapFire();
  virtual void map(float value, float& R, float& G, float& B) const;
  virtual QString getName() const;
};

inline QString ColorMapFire::getName() const
{
  return QString("Fire");
}

#endif /* FIRE_H_ */
