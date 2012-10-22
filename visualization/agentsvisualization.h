/*
 * agentsvisualization.h
 *
 *  Created on: 25-nov-2009
 *      Author: El-Kebir, M.
 */

#ifndef AGENTSVISUALIZATION_H_
#define AGENTSVISUALIZATION_H_

#include "visualization.h"
#include "gui/agentswidget.h"
#include "simulation/macrophage.h"  //Needed for Mac::NSTATES
#include <boost/multi_array.hpp>

class ScalarAgentGrid;

class AgentsVisualization : public Visualization
{
private:
  const ScalarAgentGrid* _pScalarAgentGrid;
  bool _drawGrid;
  bool _macFilter[Mac::NSTATES][AgentsWidget::NUM];
  bool _drawTgam;
  bool _drawTcyt;
  bool _drawTreg;
  bool _drawCas;
  bool _drawSrc;
  bool _drawExtMtb;
  bool _drawSquares;
  float _gridAlpha;
  float _gridHeight;
  bool _drawSrcMac;
  bool _drawSrcTgam;
  bool _drawSrcTcyt;
  bool _drawSrcTreg;
  int _drawm1m2;
  double _m1m2thres;
  int _selRow;
  int _selCol;
  boost::multi_array<float, 3> stratify;
  void drawGrid() const;
  void drawQuad(int i, int j) const;
  void drawCircle(float x, float y, float z, float r) const;
  void drawCross(int i, int j) const;
  void drawMark(int i, int j) const;
  bool drawMac(const Mac* pMac, int x, int y, int p, const double minRatio, const double maxRatio) const;
  void drawTcell(const AgentType a_t, int x, int y, int p=0) const;
  void drawCell(float i, float j, float k, float r=0.25f) const;
  void drawM1M2Agent(const Agent* a, int p, const double minRatio, const double maxRatio) const;

public:
  AgentsVisualization(int DIM, const ScalarAgentGrid* pScalarAgentGrid);
  virtual ~AgentsVisualization();
  void visualize(bool blend, const Simulation* pSimulation, const ColorMap* pColorMap) const;
  bool getDrawGrid() const;
  void setDrawGrid(bool drawGrid);
  float getGridAlpha() const;
  void setGridAlpha(float alpha);
  void setDrawTgam(bool value);
  void setDrawTcyt(bool value);
  void setDrawTreg(bool value);
  void setDrawCas(bool value);
  void setDrawSrc(bool value);
  void setDrawExtMtb(bool value);
  void setDrawSquares(bool value);
  void setGridHeight(float gridHeight);
  void setPredicates(bool mac, bool tgam, bool tcyt, bool treg);
  void setSelection(int row, int col);
  void setDrawM1M2(int value, double threshold);
  void setDrawMac(Mac::State s, AgentsWidget::MacSecondStates s2, bool enable)
  {
    _macFilter[s][s2] = enable;
  }
};

inline void AgentsVisualization::setSelection(int row, int col)
{
  if (_selRow == row && _selCol == col)
    {
      _selRow = _selCol = -1;
    }
  else
    {
      _selRow = row;
      _selCol = col;
    }
}

inline void AgentsVisualization::setPredicates(bool mac, bool tgam, bool tcyt, bool treg)
{
  _drawSrcMac = mac;
  _drawSrcTgam = tgam;
  _drawSrcTcyt = tcyt;
  _drawSrcTreg = treg;
}

inline void AgentsVisualization::setGridHeight(float gridHeight)
{
  _gridHeight = gridHeight;
}

inline void AgentsVisualization::setDrawTgam(bool value)
{
  _drawTgam = value;
}

inline void AgentsVisualization::setDrawTcyt(bool value)
{
  _drawTcyt = value;
}

inline void AgentsVisualization::setDrawTreg(bool value)
{
  _drawTreg = value;
}

inline void AgentsVisualization::setDrawCas(bool value)
{
  _drawCas = value;
}

inline void AgentsVisualization::setDrawSrc(bool value)
{
  _drawSrc = value;
}

inline void AgentsVisualization::setDrawExtMtb(bool value)
{
  _drawExtMtb = value;
}

inline void AgentsVisualization::setDrawSquares(bool value)
{
  _drawSquares = value;
}

inline void AgentsVisualization::drawCross(int i, int j) const
{
  /*glBegin(GL_LINES);
  	glVertex3f((j + .5) * deltaX, i * deltaY, _gridHeight);
  	glVertex3f((j + .5) * deltaX, (i + 1) * deltaY, _gridHeight);
  	glVertex3f((j + .1) * deltaX, (i + .65) * deltaY, _gridHeight);
  	glVertex3f((j + .9) * deltaX, (i + .65) * deltaY, _gridHeight);
  glEnd();*/
  glBegin(GL_LINES);
  glVertex3f(j * _deltaX, i * _deltaY, 0.0f);
  glVertex3f((j + 1) * _deltaX, (i + 1) * _deltaY, 0.0f);
  glVertex3f(j * _deltaX, (i + 1) * _deltaY, 0.0f);
  glVertex3f((j + 1) * _deltaX, i * _deltaY, 0.0f);
  glEnd();

}

inline void AgentsVisualization::drawMark(int i, int j) const
{
  glBegin(GL_LINES);
  glVertex3f((j + 0.5) * _deltaX, i * _deltaY, 0.01f);
  glVertex3f((j + 0.5) * _deltaX, (i + 1) * _deltaY, 0.01f);
  glVertex3f(j * _deltaX, (i + 0.5) * _deltaY, 0.01f);
  glVertex3f((j + 1) * _deltaX, (i + 0.5) * _deltaY, 0.01f);
  glEnd();
}

inline void AgentsVisualization::drawQuad(int i, int j) const
{
  glBegin(GL_QUADS);
  glVertex3f(j * _deltaX, i * _deltaY, _gridHeight);
  glVertex3f(j * _deltaX, (i + 1) * _deltaY, _gridHeight);
  glVertex3f((j + 1) * _deltaX, (i + 1) * _deltaY, _gridHeight);
  glVertex3f((j + 1) * _deltaX, i * _deltaY, _gridHeight);
  glEnd();
}
inline void AgentsVisualization::drawCircle(float x, float y, float z, float r) const
{
  glBegin(GL_POLYGON);
    float angle = 0.0f;
    for(int i=0;i<15; i++) {
      angle = i*2*M_PI/15.0f;
      glVertex3f(x + cos(angle) * r, y + sin(angle) * r, z);
    }
  glEnd();
}

inline void AgentsVisualization::drawCell(float i, float j, float k, float R) const
{
  assert(R<=0.5f);  //Any larger and the cell overlaps the adjacent square
  static float currentColor[4];
  glGetFloatv(GL_CURRENT_COLOR, currentColor);
  drawCircle((j+0.5f)*_deltaX, (i+0.5f)*_deltaY, k+_gridHeight+0.001f, R*_deltaX);
  for(int z=0;z<3;z++)
    currentColor[z] *= 0.25f;  //Make 75% darker for nucleus
  glColor4fv(currentColor);
  float rngx = stratify[floor(j)][floor(i)][0];
  float rngy = stratify[floor(j)][floor(i)][1];
  drawCircle((j+(rngx*R+1)*0.5f)*_deltaX, (i+(rngy*R+1)*0.5f)*_deltaY, k+_gridHeight+0.002f, R*0.5f*_deltaX);
}

inline bool AgentsVisualization::getDrawGrid() const
{
  return _drawGrid;
}

inline void AgentsVisualization::setDrawGrid(bool drawGrid)
{
  _drawGrid = drawGrid;
}

inline float AgentsVisualization::getGridAlpha() const
{
  return _gridAlpha;
}

inline void AgentsVisualization::setGridAlpha(float alpha)
{
  _gridAlpha = alpha;
}
inline void AgentsVisualization::setDrawM1M2(int value, double threshold)
{
  _drawm1m2 = value;
  _m1m2thres = threshold;
}


#endif /* AGENTSVISUALIZATION_H_ */
