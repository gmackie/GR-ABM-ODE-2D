/*
 * glyphcone.h
 *
 *  Created on: 22-sep-2008
 *      Author: S030858
 */

#ifndef GLYPHCONE_H_
#define GLYPHCONE_H_

#include "glyph.h"
#include "vectordatasets/vector.h"

class GlyphCone: public Glyph
{
private:
	GLUquadric* _pConeQuadric;
	GLUquadric* _pConeDiscQuadric;

public:
	GlyphCone();
	virtual ~GlyphCone();
	virtual void draw(const vec2f& origin, const vec2f& vec) const;
};

inline GlyphCone::GlyphCone()
	: Glyph()
	, _pConeQuadric(gluNewQuadric())
	, _pConeDiscQuadric(gluNewQuadric())
{
}

inline GlyphCone::~GlyphCone()
{
	gluDeleteQuadric(_pConeQuadric);
	gluDeleteQuadric(_pConeDiscQuadric);
}

inline void GlyphCone::draw(const vec2f& origin, const vec2f& vec) const
{
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    glPushMatrix();
    float magVec = vec.magnitude();
    setTransAndRotateMatrix(origin, vec);
    gluCylinder(_pConeQuadric, magVec*0.3, 0.0, magVec, 5, 1);
    gluDisk(_pConeDiscQuadric,0,magVec*0.3,5,1);
	glPopMatrix();

	glDisable(GL_COLOR_MATERIAL);
	glDisable(GL_LIGHTING);
    glDisable(GL_LIGHT0);
}

#endif /* GLYPHCONE_H_ */
