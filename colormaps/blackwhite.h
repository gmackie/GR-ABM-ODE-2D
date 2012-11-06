/*
 * blackwhite.h
 *
 *  Created on: 8-sep-2008
 *      Author: s030858
 */

#ifndef BLACKWHITE_H_
#define BLACKWHITE_H_

#include "colormap.h"
#include <QString>

class ColorMapBlackWhite : public ColorMap
{
public:
  ColorMapBlackWhite();
  ~ColorMapBlackWhite();
  void map(float value, float& R, float& G, float& B) const;
  QString getName() const;
  COLORMAP getType() const;
};

inline QString ColorMapBlackWhite::getName() const
{
  return QString("Grayscale");
}

inline COLORMAP ColorMapBlackWhite::getType() const
{
  return CMAPGRAY;
}

#endif /* BLACKWHITE_H_ */
