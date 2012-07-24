/*
 * isolinesvisualization.cpp
 *
 *  Created on: 7-nov-2008
 *      Author: s030858
 */

#include <assert.h>
#include "isolinesvisualization.h"
#include "scalardatasets/scalarnormalizer.h"
#include "scalardatasets/scalargrid.h"

/* Marching triangles
 *
 * C
 * | \
 * |  \
 * |   \
 * A --- B
 *
 * 0 : <  target value
 * 1 : >= target value
 *
 * A B C	AB BC AC
 * =====	========
 * 0 0 0	-  -  -
 * 0 0 1	-  *  *
 * 0 1 0	*  *  -
 * 0 1 1	*  -  *
 * 1 0 0	*  -  *
 * 1 0 1	*  *  -
 * 1 1 0	-  *  *
 * 1 1 1	-  -  -
 */

IsolinesVisualization::IsolinesVisualization(int DIM, const ScalarNormalizer* pScalarNormalizer, const ScalarGrid* pScalarGrid)
  : Visualization(DIM)
  , _pScalarNormalizer(pScalarNormalizer)
  , _pScalarGrid(pScalarGrid)
  , _targetValueVector(1, _ISOLINES_TARGET_VALUE)
  , _lineWidth(_ISOLINES_LINE_WIDTH)
{
  assert(pScalarNormalizer && pScalarGrid);
}

IsolinesVisualization::~IsolinesVisualization()
{
}

void IsolinesVisualization::processTriangle(const vec2f& a, const vec2f& b, const vec2f& c,
    float val_a, float val_b, float val_c, float targetValue) const
{
  int bitarray = 0;

  if (val_a >= targetValue)
    bitarray |= 4;
  if (val_b >= targetValue)
    bitarray |= 2;
  if (val_c >= targetValue)
    bitarray |= 1;

  /*
  std::cout << "Bitarray: " << bitarray
  	<< "\ttarget: " << _targetValue
  	<< "\ta: " << val_a
  	<< "\tb: " << val_b
  	<< "\tc: " << val_c
  	<< std::endl;
  */

  if (bitarray != 0 && bitarray != 7)
    {
      vec2f v0, v1;
      if (bitarray == 1 || bitarray == 6)
        {
          v0 = interpolate(b, c, val_b, val_c, targetValue);
          v1 = interpolate(a, c, val_a, val_c, targetValue);
        }
      else if (bitarray == 2 || bitarray == 5)
        {
          v0 = interpolate(a, b, val_a, val_b, targetValue);
          v1 = interpolate(b, c, val_b, val_c, targetValue);
        }
      else if (bitarray == 3 || bitarray == 4)
        {
          v0 = interpolate(a, b, val_a, val_b, targetValue);
          v1 = interpolate(a, c, val_a, val_c, targetValue);
        }

      /*
      std::cout << "v0 = (" << v0[0] << "," << v0[1] << ")"
      	<< "\tv1 = (" << v1[0] << "," << v1[1] << ")"
      	<< std::endl;
      */

      v0[0] *= _deltaX;
      v0[1] *= _deltaY;
      v1[0] *= _deltaX;
      v1[1] *= _deltaY;

      glBegin(GL_LINES);
      glVertex2fv(v0.v());
      glVertex2fv(v1.v());
      glEnd();
    }
}

vec2f IsolinesVisualization::interpolate(const vec2f& p_i, const vec2f& p_j, float v_i, float v_j, float v) const
{
  float t = (v - v_i) / (v_j - v_i);
  vec2f res = p_i + t * (p_j - p_i);

  /*
  std::cout << "p_i = (" << p_i[0] << "," << p_i[1] << ")"
  	<< "\tp_j = (" << p_j[0] << "," << p_j[1] << ")"
  	<< "\tt: " << t
  	<< "\tres = (" << res[0] << "," << res[1] << ")" << std::endl;
   */
  return res;
}

void IsolinesVisualization::visualize(bool blend, const Simulation*, const ColorMap* pColorMap) const
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

  const std::vector<ScalarGridItem>& grid = _pScalarGrid->getGrid();

  assert(pColorMap);
  glLineWidth(_lineWidth);

  for (int i = 0; i < _DIM - 1; i++)
    {
      for (int j = 0; j < _DIM - 1; j++)
        {
          // c d
          // a b
          const ScalarGridItem& gridItemA = grid[i * _DIM + j];
          const ScalarGridItem& gridItemB = grid[i * _DIM + j + 1];
          const ScalarGridItem& gridItemC = grid[(i + 1) * _DIM + j];
          const ScalarGridItem& gridItemD = grid[(i + 1) * _DIM + j + 1];

          vec2f a = gridItemA.pos;
          vec2f b = gridItemB.pos;
          vec2f c = gridItemC.pos;
          vec2f d = gridItemD.pos;

          for (size_t k = 0; k < _targetValueVector.size(); k++)
            {
              float normalizedValue = _targetValueVector[k];
              float value = _pScalarNormalizer->denormalize(normalizedValue);
              applyColorMap(pColorMap, normalizedValue);

              processTriangle(a, b, c, gridItemA.scalar, gridItemB.scalar,
                              gridItemC.scalar, value);
              processTriangle(b, c, d, gridItemB.scalar, gridItemC.scalar,
                              gridItemD.scalar, value);
            }
        }
    }

  glLineWidth(1.0f);
  glDisable(GL_BLEND);
  glEnable(GL_DEPTH_TEST);
}
