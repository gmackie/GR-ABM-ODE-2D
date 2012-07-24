/*
 * heightplotvisualization.cpp
 *
 *  Created on: 18-nov-2008
 *      Author: s030858
 */

#include "heightplotvisualization.h"
#include "scalardatasets/scalargrid.h"
#include "scalardatasets/scalarnormalizer.h"
#include "vectordatasets/vectorgradientdataset.h"

HeightPlotVisualization::HeightPlotVisualization(int DIM, const ScalarNormalizer* pScalarHeightNormalizer,
    const ScalarNormalizer* pScalarHeightColorNormalizer,
    const ScalarGrid* pScalarHeightGrid, const ScalarGrid* pScalarHeightColorGrid,
    VectorGradientDataset* pGradientDataset)
  : Visualization(DIM)
  , _pScalarHeightNormalizer(pScalarHeightNormalizer)
  , _pScalarHeightColorNormalizer(pScalarHeightColorNormalizer)
  , _pScalarHeightGrid(pScalarHeightGrid)
  , _pScalarHeightColorGrid(pScalarHeightColorGrid)
  , _maxHeight(_HEIGHTPLOT_MAX_HEIGHT)
  , _drawGrid(true)
  , _gridHeight(_HEIGHTPLOT_GRID_HEIGHT)
  , _gridAlpha(_HEIGHTPLOT_GRID_ALPHA)
  , _pGradientDataset(pGradientDataset)
{
  assert(pScalarHeightNormalizer && pScalarHeightNormalizer &&
         pScalarHeightGrid && pScalarHeightColorGrid && pGradientDataset);
}

HeightPlotVisualization::~HeightPlotVisualization()
{
  delete _pGradientDataset;
}

void HeightPlotVisualization::visualize(bool blend, const Simulation* pSimulation, const ColorMap* pColorMap) const
{
  if (blend)
    {
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE);
      //glDisable(GL_DEPTH_TEST);
    }
  else
    {
      glDisable(GL_BLEND);
      //glEnable(GL_DEPTH_TEST);
    }

  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);

  const std::vector<ScalarGridItem>& heightGrid = _pScalarHeightGrid->getGrid();
  const std::vector<ScalarGridItem>& heightColorGrid = _pScalarHeightColorGrid->getGrid();

  assert(pColorMap);

  glPolygonMode(GL_FRONT, GL_FILL);

  if (_drawGrid)
    {
      drawGrid();
    }

  vec2f tempVec;

  glBegin(GL_TRIANGLES);
  for (int i = 0; i < _DIM - 1; i++)
    {
      for (int j = 0; j < _DIM - 1; j++)
        {
          // c d
          // a b
          const ScalarGridItem& heightGridItemA = heightGrid[i * _DIM + j];
          const ScalarGridItem& heightGridItemB = heightGrid[i * _DIM + j + 1];
          const ScalarGridItem& heightGridItemC = heightGrid[(i + 1) * _DIM + j];
          const ScalarGridItem& heightGridItemD = heightGrid[(i + 1) * _DIM + j + 1];

          const ScalarGridItem& heightColorGridItemA = heightColorGrid[i * _DIM + j];
          const ScalarGridItem& heightColorGridItemB = heightColorGrid[i * _DIM + j + 1];
          const ScalarGridItem& heightColorGridItemC = heightColorGrid[(i + 1) * _DIM + j];
          const ScalarGridItem& heightColorGridItemD = heightColorGrid[(i + 1) * _DIM + j + 1];

          vec3f a, normalA;
          a[0] = heightGridItemA.pos[0];
          a[1] = heightGridItemA.pos[1];
          a[2] = _pScalarHeightNormalizer->normalize(heightGridItemA.scalar) * _maxHeight;

          _pGradientDataset->getVectorBL(pSimulation, heightColorGridItemA.pos, tempVec);
          normalA[0] = tempVec[0];
          normalA[1] = tempVec[1];
          normalA[2] = 1.0f;
          normalA.normalize();

          vec3f b, normalB;
          b[0] = heightGridItemB.pos[0];
          b[1] = heightGridItemB.pos[1];
          b[2] = _pScalarHeightNormalizer->normalize(heightGridItemB.scalar) * _maxHeight;

          _pGradientDataset->getVectorBL(pSimulation, heightColorGridItemB.pos, tempVec);
          normalB[0] = tempVec[0];
          normalB[1] = tempVec[1];
          normalB[2] = 1.0f;
          normalB.normalize();

          vec3f c, normalC;
          c[0] = heightGridItemC.pos[0];
          c[1] = heightGridItemC.pos[1];
          c[2] = _pScalarHeightNormalizer->normalize(heightGridItemC.scalar) * _maxHeight;

          _pGradientDataset->getVectorBL(pSimulation, heightColorGridItemC.pos, tempVec);
          normalC[0] = tempVec[0];
          normalC[1] = tempVec[1];
          normalC[2] = 1.0f;
          normalC.normalize();

          vec3f d, normalD;
          d[0] = heightGridItemD.pos[0];
          d[1] = heightGridItemD.pos[1];
          d[2] = _pScalarHeightNormalizer->normalize(heightGridItemD.scalar) * _maxHeight;

          _pGradientDataset->getVectorBL(pSimulation, heightColorGridItemD.pos, tempVec);
          normalD[0] = tempVec[0];
          normalD[1] = tempVec[1];
          normalD[2] = 1.0f;
          normalD.normalize();

          a[0] *= _deltaX;
          a[1] *= _deltaY;
          b[0] *= _deltaX;
          b[1] *= _deltaY;
          c[0] *= _deltaX;
          c[1] *= _deltaY;
          d[0] *= _deltaX;
          d[1] *= _deltaY;

          glNormal3fv(normalA.v());
          applyColorMap(pColorMap, _pScalarHeightColorNormalizer->normalize(heightColorGridItemA.scalar));
          glVertex3fv(a.v());

          glNormal3fv(normalB.v());
          applyColorMap(pColorMap, _pScalarHeightColorNormalizer->normalize(heightColorGridItemB.scalar));
          glVertex3fv(b.v());

          glNormal3fv(normalC.v());
          applyColorMap(pColorMap, _pScalarHeightColorNormalizer->normalize(heightColorGridItemC.scalar));
          glVertex3fv(c.v());

          glNormal3fv(normalB.v());
          applyColorMap(pColorMap, _pScalarHeightColorNormalizer->normalize(heightColorGridItemB.scalar));
          glVertex3fv(b.v());

          glNormal3fv(normalD.v());
          applyColorMap(pColorMap, _pScalarHeightColorNormalizer->normalize(heightColorGridItemD.scalar));
          glVertex3fv(d.v());

          glNormal3fv(normalC.v());
          applyColorMap(pColorMap, _pScalarHeightColorNormalizer->normalize(heightColorGridItemC.scalar));
          glVertex3fv(c.v());
        }
    }
  glEnd();

  glDisable(GL_COLOR_MATERIAL);
  glDisable(GL_LIGHTING);
  glDisable(GL_LIGHT0);
  glDisable(GL_BLEND);
}

void HeightPlotVisualization::drawGrid() const
{
  int nrSamplesX = std::min(_DIM, 75);
  int nrSamplesY = std::min(_DIM, 75);

  float deltaX = _MAX_X / (nrSamplesX - 1);
  float deltaY = _MAX_Y / (nrSamplesY - 1);

  glColor4f(1.0f, 1.0f, 1.0f, _gridAlpha);

  glBegin(GL_LINES);
  for (int i = 0; i < nrSamplesX; i++)
    {
      glVertex3f(i * deltaX, _MIN_Y, _gridHeight);
      glVertex3f(i * deltaX, _MAX_Y, _gridHeight);
    }
  for (int i = 0; i < nrSamplesY; i++)
    {
      glVertex3f(_MIN_X, i * deltaY, _gridHeight);
      glVertex3f(_MAX_X, i * deltaY, _gridHeight);
    }
  glEnd();
}

void HeightPlotVisualization::setGradientDataset(VectorGradientDataset* pGradientDataset)
{
  delete _pGradientDataset;
  _pGradientDataset = pGradientDataset;
}
