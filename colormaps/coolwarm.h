/*
 * coolwarm.h
 *
 *  Created on: 12-sep-2008
 *      Author: s030858
 */

#ifndef COOLWARM_H_
#define COOLWARM_H_

#include "colormap.h"
#include <QString>

class ColorMapCoolWarm : public ColorMap
{
public:
  ColorMapCoolWarm();
  virtual ~ColorMapCoolWarm();
  virtual void map(float value, float& R, float& G, float& B) const;
  virtual QString getName() const;
};

inline QString ColorMapCoolWarm::getName() const
{
  return QString("Cool/warm");
}

#endif /* COOLWARM_H_ */
