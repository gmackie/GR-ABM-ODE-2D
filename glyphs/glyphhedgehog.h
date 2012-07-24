/*
 * glyphhedgehog.h
 *
 *  Created on: 22-sep-2008
 *      Author: S030858
 */

#ifndef GLYPHHEDGEHOG_H_
#define GLYPHHEDGEHOG_H_

#include "glyph.h"
#include "vectordatasets/vector.h"

class GlyphHedgehog : public Glyph
{
public:
  GlyphHedgehog();
  virtual ~GlyphHedgehog();
  virtual void draw(const vec2f& origin, const vec2f& vec) const;
};


inline GlyphHedgehog::GlyphHedgehog()
  : Glyph()
{
}

inline GlyphHedgehog::~GlyphHedgehog()
{
}

inline void GlyphHedgehog::draw(const vec2f& origin, const vec2f& vec) const
{
  glBegin(GL_LINES);
  glVertex2fv(origin.v());
  glVertex2fv((origin + vec).v());
  glEnd();
}

#endif /* GLYPHHEDGEHOG_H_ */
