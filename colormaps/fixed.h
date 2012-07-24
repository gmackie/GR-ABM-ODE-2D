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
  virtual ~ColorMapFixed();
  virtual void map(float value, float& R, float& G, float& B) const;
  virtual QString getName() const;
};

inline QString ColorMapFixed::getName() const
{
  return QString("Fixed");
}

#endif /* FIXED_H_ */
