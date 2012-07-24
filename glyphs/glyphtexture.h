/*
 * glyphtexture.h
 *
 *  Created on: 6-okt-2008
 *      Author: s020751
 */

#ifndef GLYPHTEXTURE_H_
#define GLYPHTEXTURE_H_

#include "glyph.h"
#include "vectordatasets/vector.h"
#include <QImage>
#include <QString>
#include <QGLWidget>

class GlyphTexture : public Glyph
{
private:
  QImage _texture;
  GLuint _texName;
  float _deltaX;
  float _deltaY;

public:
  GlyphTexture(const QString& fileName);
  virtual ~GlyphTexture();
  virtual void draw(const vec2f& origin, const vec2f& vec) const;
  void setDeltaY(float deltaY);
};

inline void GlyphTexture::setDeltaY(float deltaY)
{
  _deltaY = deltaY;
}

inline GlyphTexture::GlyphTexture(const QString& fileName)
  : Glyph()
{
  _texture = QGLWidget::convertToGLFormat( QImage(fileName) );

  glGenTextures(1, &_texName);
  glBindTexture(GL_TEXTURE_2D, _texName);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);

  glTexImage2D(GL_TEXTURE_2D, 0, 4, _texture.width(), _texture.height(), 0, GL_RGBA,
               GL_UNSIGNED_BYTE, _texture.bits());
}

inline GlyphTexture::~GlyphTexture()
{
}

inline void GlyphTexture::draw(const vec2f& origin, const vec2f& vec) const
{
  glDisable(GL_COLOR_MATERIAL);
  glEnable(GL_TEXTURE_2D);
  glAlphaFunc(GL_GREATER, 0.1f);
  glEnable(GL_ALPHA_TEST);

  glPushMatrix();
  setTransAndRotateMatrix(origin, vec);

  float mag = vec.magnitude();
  float a = 0.5f * sqrt(2.0) * mag;
  float y = sqrt(a * a - 0.25f * mag * mag);

  glBegin(GL_QUADS);
  glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
  glTexCoord2d(0.0, 0.0);
  glVertex3f(0.0f, 0.0f, 0.0f);
  glTexCoord2d(1.0, 0.0);
  glVertex3f(0.0f, -y, 0.5f * mag);
  glTexCoord2d(1.0, 1.0);
  glVertex3f(0.0f, 0.0f, mag);
  glTexCoord2d(0.0, 1.0);
  glVertex3f(0.0f, y, 0.5f * mag);
  glEnd();

  glPopMatrix();
  glDisable(GL_ALPHA_TEST);
  glDisable(GL_TEXTURE_2D);
}

#endif /* GLYPHTEXTURE_H_ */

