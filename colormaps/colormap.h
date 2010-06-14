/*
 * colormapping.h
 *
 *  Created on: 8-sep-2008
 *      Author: s030858
 */

#ifndef COLORMAPPING_H_
#define COLORMAPPING_H_

#include <QString>
#include "grviz.h"

class ColorMap
{
protected:
	bool _invert;
	int _nrBands;
	float _hueDelta;
	float _satDelta;
	float _valDelta;
	float _alpha;

	void bin(float& value) const;
	void invert(float& value) const;
	void applyDeltas(float& r, float& g, float& b) const;
	void rgb2hsv(float r, float g, float b, float& h, float& s, float& v) const;
	void hsv2rgb(float h, float s, float v, float& r, float& g, float& b) const;

public:
	ColorMap();
	virtual ~ColorMap();
	virtual void map(float value, float& R, float& G, float& B) const = 0;
	int getNrBands() const;
	void setNrBands(int nrBands);
	void setHueDelta(float hueDelta);
	float getHueDelta() const;
	void setSatDelta(float satDelta);
	float getSatDelta() const;
	void setValDelta(float valDelta);
	float getValDelta() const;
	void setAlpha(float alpha);
	float getAlpha() const;
	void setInvert(bool invert);
	bool getInvert() const;
	virtual QString getName() const = 0;
};

inline void ColorMap::setInvert(bool invert)
{
	_invert = invert;
}

inline bool ColorMap::getInvert() const
{
	return _invert;
}

inline void ColorMap::invert(float& value) const
{
	value = 1 - value;
}

inline int ColorMap::getNrBands() const
{
	return _nrBands;
}

inline void ColorMap::setNrBands(int nrBands)
{
	_nrBands = nrBands;
}

inline void ColorMap::setHueDelta(float hueDelta)
{
	_hueDelta = hueDelta;
}

inline float ColorMap::getHueDelta() const
{
	return _hueDelta;
}

inline void ColorMap::setSatDelta(float satDelta)
{
	_satDelta = satDelta;
}

inline float ColorMap::getSatDelta() const
{
	return _satDelta;
}

inline void ColorMap::setValDelta(float valDelta)
{
	_valDelta = valDelta;
}

inline float ColorMap::getValDelta() const
{
	return _valDelta;
}

inline void ColorMap::setAlpha(float alpha)
{
	_alpha = alpha;
}

inline float ColorMap::getAlpha() const
{
	return _alpha;
}

#endif /* COLORMAPPING_H_ */
