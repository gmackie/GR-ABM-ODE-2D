/*
 * glyphvisualization.h
 *
 *  Created on: 21-sep-2008
 *      Author: S030858
 */

#ifndef GLYPHVISUALIZATION_H_
#define GLYPHVISUALIZATION_H_

#include "visualization.h"

class ScalarNormalizer;
class ScalarGrid;
class VectorGrid;
class Glyph;

class GlyphVisualization : public Visualization
{
public:
  GlyphVisualization(int DIM, const ScalarNormalizer* pScalarNormalizer,
                     const ScalarGrid* pScalarGrid, const VectorGrid* pVectorGrid);
  virtual ~GlyphVisualization();
  virtual void visualize(bool blend, const Simulation* pSimulation, const ColorMap* pColorMap) const;
  void setScaleFactor(int value);
  void setGlyph(Glyph* pGlyph);
  Glyph* getGlyph() const;
  void setClamping(bool clamp);

private:
  const ScalarNormalizer* _pScalarNormalizer;
  const ScalarGrid* _pScalarGrid;
  const VectorGrid* _pVectorGrid;
  int _scaleFactor;
  Glyph* _pGlyph;
  bool _clamp;
};

inline Glyph* GlyphVisualization::getGlyph() const
{
  return _pGlyph;
}

inline void GlyphVisualization::setScaleFactor(int value)
{
  _scaleFactor = value;
}

inline void GlyphVisualization::setClamping(bool clamp)
{
  _clamp = clamp;
}

#endif /* GLYPHVISUALIZATION_H_ */
