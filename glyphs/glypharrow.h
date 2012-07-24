/*
 * glypharrow.h
 *
 *  Created on: 25-sep-2008
 *      Author: S030858
 */

#ifndef GLYPHARROW_H_
#define GLYPHARROW_H_

#include "glyph.h"

class GlyphArrow: public Glyph
{
public:
  GlyphArrow();
  virtual ~GlyphArrow();
  virtual void draw(const vec2f& origin, const vec2f& vec) const;

private:
  GLUquadric* _pCylinderQuadric;
  GLUquadric* _pCylinderDiscQuadric;
  GLUquadric* _pConeQuadric;
  GLUquadric* _pConeDiscQuadric;
};

inline GlyphArrow::GlyphArrow()
  : Glyph()
  , _pCylinderQuadric(gluNewQuadric())
  , _pCylinderDiscQuadric(gluNewQuadric())
  , _pConeQuadric(gluNewQuadric())
  , _pConeDiscQuadric(gluNewQuadric())
{
}

inline GlyphArrow::~GlyphArrow()
{
  gluDeleteQuadric(_pCylinderQuadric);
  gluDeleteQuadric(_pCylinderDiscQuadric);
  gluDeleteQuadric(_pConeQuadric);
  gluDeleteQuadric(_pConeDiscQuadric);
}

inline void GlyphArrow::draw(const vec2f& origin, const vec2f& vec) const
{
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);

  float magVec = vec.magnitude();

  float cylinderDiameter = magVec * 0.2f;
  float cone = magVec * 0.4f;

  glPushMatrix();
  setTransAndRotateMatrix(origin, vec);
  gluCylinder(_pCylinderQuadric, cylinderDiameter, cylinderDiameter, magVec - cone, 5, 1);
  gluDisk(_pCylinderDiscQuadric,0,cylinderDiameter,5,1);
  glPopMatrix();

  vec2f newVec = vec;
  newVec *= (magVec - cone) / magVec;

  glPushMatrix();
  setTransAndRotateMatrix(origin + newVec, newVec);
  gluCylinder(_pConeQuadric, cone, 0.0, cone, 5, 1);
  gluDisk(_pConeDiscQuadric, 0, cone, 5, 1);
  glPopMatrix();

  glDisable(GL_COLOR_MATERIAL);
  glDisable(GL_LIGHTING);
  glDisable(GL_LIGHT0);
}

#endif /* GLYPHARROW_H_ */
