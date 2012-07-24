/*
 * smokevisualization.cpp
 *
 *  Created on: 21-sep-2008
 *      Author: S030858
 */

#include "smokevisualization.h"
#include "scalardatasets/scalarnormalizer.h"
#include "scalardatasets/scalargrid.h"

SmokeVisualization::SmokeVisualization(int DIM, const ScalarNormalizer* pScalarNormalizer, const ScalarGrid* pScalarGrid)
  : Visualization(DIM)
  , _pScalarNormalizer(pScalarNormalizer)
  , _pScalarGrid(pScalarGrid)
{
  assert(pScalarNormalizer && pScalarGrid);
}

SmokeVisualization::~SmokeVisualization()
{
}

void SmokeVisualization::visualize(bool blend, const Simulation*, const ColorMap* pColorMap) const
{
  if (blend)
    {
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    }

  const std::vector<ScalarGridItem>& grid = _pScalarGrid->getGrid();

  assert(pColorMap);

  double px, py;
  float value;

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  for (int i = 0; i < _DIM - 1; i++) //draw smoke
    {
      glBegin(GL_TRIANGLE_STRIP);

      // . .
      // * .
      ScalarGridItem gridItem = grid[i * _DIM];
      px = gridItem.pos[0] * _deltaX;
      py = gridItem.pos[1] * _deltaY;

      value = _pScalarNormalizer->normalize(gridItem.scalar);

      applyColorMap(pColorMap, value);
      glVertex2f(px, py);

      for (int j = 0; j < _DIM - 1; j++)
        {
          // * .
          // . .
          gridItem = grid[j + (i+1) * _DIM];

          value = _pScalarNormalizer->normalize(gridItem.scalar);
          px = gridItem.pos[0] * _deltaX;
          py = gridItem.pos[1] * _deltaY;
          applyColorMap(pColorMap, value);
          glVertex2f(px, py);

          // . .
          // . *
          gridItem = grid[j + 1 + (i) * _DIM];

          value = _pScalarNormalizer->normalize(gridItem.scalar);
          px = gridItem.pos[0] * _deltaX;
          py = gridItem.pos[1] * _deltaY;
          applyColorMap(pColorMap, value);
          glVertex2f(px, py);
        }

      // . *
      // . .
      gridItem = grid[_DIM - 1 + (i + 1) * _DIM];

      value = _pScalarNormalizer->normalize(gridItem.scalar);
      px = gridItem.pos[0] * _deltaX;
      py = gridItem.pos[1] * _deltaY;
      applyColorMap(pColorMap, value);
      glVertex2f(px, py);
      glEnd();
    }

  glDisable(GL_BLEND);
}
