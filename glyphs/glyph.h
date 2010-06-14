/*
 * glyph.h
 *
 *  Created on: 22-sep-2008
 *      Author: S030858
 */

#ifndef GLYPH_H_
#define GLYPH_H_

#include "grviz.h"
#include "colormaps/colormap.h"
#include "vectordatasets/vector.h"

class Glyph
{
public:
	Glyph();
	virtual ~Glyph();
	virtual void draw(const vec2f& origin, const vec2f& vec) const = 0;

protected:
	virtual void setTransAndRotateMatrix(const vec2f& origin, const vec2f& vec) const;
};

inline Glyph::Glyph()
{
}

inline Glyph::~Glyph()
{
}

inline void Glyph::setTransAndRotateMatrix(const vec2f& origin, const vec2f& vec) const
{
	//glTranslatef(origin[0], origin[1], 0.0);
	//glRotatef(90.0, 0.0, 1.0, 0.0);
	//glRotatef(-45.0, 1.0, 0.0, 0.0);

    float magVec = vec.magnitude();

    // matrices in opengl are in column major order
    float m[16];
    m[0] = m[1] = m[3] = m[6] = m[7] = m[10] = m[11] = m[14] = 0;
    m[2] = -1;
    m[4] = -1 * vec[1] / magVec; // sin(-a) = -sin(a)
    m[5] = m[8] = vec[0] / magVec; // cos(-a) = cos(a)
    m[9] = -m[4];
    m[12] = origin[0];
    m[13] = origin[1];
    m[15] = 1;

    glMultMatrixf(m);
}

#endif /* GLYPH_H_ */
