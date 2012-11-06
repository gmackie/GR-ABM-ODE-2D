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
  ~ColorMapCoolWarm();
  void map(float value, float& R, float& G, float& B) const;
  QString getName() const;
  COLORMAP getType() const;
};

inline QString ColorMapCoolWarm::getName() const
{
  return QString("Cool/warm");
}

inline COLORMAP ColorMapCoolWarm::getType() const
{
  return CMAPCOOLWARM;
}

#endif /* COOLWARM_H_ */
