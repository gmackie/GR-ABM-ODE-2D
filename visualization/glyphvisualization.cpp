/*
 * glyphvisualization.cpp
 *
 *  Created on: 21-sep-2008
 *      Author: S030858
 */

#include "glyphvisualization.h"
#include <vector>
#include <iostream>
#include "vectordatasets/vectorgrid.h"
#include "scalardatasets/scalargrid.h"
#include "scalardatasets/scalarnormalizer.h"
#include "glyphs/glyphhedgehog.h"
#include "glyphs/glyphcone.h"
#include "glyphs/glyphtexture.h"
#include <iostream>
#include <math.h>

GlyphVisualization::GlyphVisualization(int DIM, const ScalarNormalizer* pScalarNormalizer,
                                       const ScalarGrid* pScalarGrid, const VectorGrid* pVectorGrid)
  : Visualization(DIM)
  , _pScalarNormalizer(pScalarNormalizer)
  , _pScalarGrid(pScalarGrid)
  , _pVectorGrid(pVectorGrid)
  , _scaleFactor(1)
  , _pGlyph(new GlyphHedgehog())
  , _clamp(true)
{
  _deltaX = _MAX_X / (_DIM);
  _deltaY = _MAX_Y / (_DIM);
}

GlyphVisualization::~GlyphVisualization()
{
  delete _pGlyph;
}

void GlyphVisualization::visualize(bool blend, const Simulation*, const ColorMap* pColorMap) const
{
  if (blend)
    {
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE);
      glDisable(GL_DEPTH_TEST);
    }
  else
    {
      glDisable(GL_BLEND);
      glEnable(GL_DEPTH_TEST);
    }

  assert(pColorMap);

  const std::vector<VectorGridItem>& vectorGrid = _pVectorGrid->getGrid();
  const std::vector<ScalarGridItem>& scalarGrid = _pScalarGrid->getGrid();

  int n = (int)vectorGrid.size();
  for (int i = 0; i < n; i++)
    {
      GlyphTexture* pGlyphTexture = dynamic_cast<GlyphTexture*>(_pGlyph);
      if (pGlyphTexture == NULL)
        {
          float scalar = _pScalarNormalizer->normalize(scalarGrid[i].scalar);
          applyColorMap(pColorMap, scalar);
        }
      else
        {
          pGlyphTexture->setDeltaY(_deltaY);
        }
      const VectorGridItem& item = vectorGrid[i];

      vec2f origin = item.origin;
      origin[0] *= _deltaX;
      origin[0] += .5f * _deltaX;
      origin[1] *= _deltaY;
      origin[1] += .5f * _deltaX;
      vec2f dest = item.vector;
      dest *= _scaleFactor / 2.0f;

      // clamp
      float radius = std::min(_deltaX, _deltaY);
      if (_clamp && dest.magnitude() > radius)
        {
          dest *= (radius / dest.magnitude());
        }

      _pGlyph->draw(origin, dest);
    }

  glDisable(GL_BLEND);
  glEnable(GL_DEPTH_TEST);
}

void GlyphVisualization::setGlyph(Glyph* pGlyph)
{
  delete _pGlyph;
  _pGlyph = pGlyph;
}
