/*
 * fixed.h
 *
 *  Created on: Jan 19, 2010
 *      Author: mohammed
 */

#ifndef FIXED_H_
#define FIXED_H_

#include "colormap.h"
#include <QColor>

class ColorMapFixed : public ColorMap
{
private:
  QColor _color;

public:
  ColorMapFixed(const QColor& color);
  ~ColorMapFixed();
  void map(float value, float& R, float& G, float& B) const;
  QString getName() const;
  COLORMAP getType() const;
};

inline QString ColorMapFixed::getName() const
{
  return QString("Fixed");
}

inline COLORMAP ColorMapFixed::getType() const
{
  return CMAPFIXED;
}

#endif /* FIXED_H_ */
