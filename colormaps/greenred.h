/*
 * greenred.h
 *
 *  Created on: 17-dec-2008
 *      Author: S020751
 */

#ifndef GREENRED_H_
#define GREENRED_H_

#include "colormap.h"
#include <QString>

class ColorMapGreenRed : public ColorMap
{
public:
	ColorMapGreenRed();
	virtual ~ColorMapGreenRed();
	virtual void map(float value, float& R, float& G, float& B) const;
	virtual QString getName() const;
};

inline QString ColorMapGreenRed::getName() const
{
	return QString("Green/Red");
}
#endif /* GREENRED_H_ */
