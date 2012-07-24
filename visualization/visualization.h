/*
 * visualization.h
 *
 *  Created on: 3-sep-2008
 *      Author: S030858
 */

#ifndef VISUALIZATION_H_
#define VISUALIZATION_H_

#include "grviz.h"
#include "colormaps/colormap.h"
#include "vectordatasets/vector.h"

class Simulation;

class Visualization
{
public:
  Visualization(int DIM);
  virtual ~Visualization();
  virtual void visualize (bool blend, const Simulation* pSimulation, const ColorMap* pColorMap) const = 0;

protected:
  virtual void applyColorMap(const ColorMap* pColorMap, float value) const;
  vec3f crossproduct(const vec3f& vec1, const vec3f& vec2) const;
  int _DIM;
  float _deltaX;
  float _deltaY;

};

inline Visualization::Visualization(int DIM)
  : _DIM(DIM)
  , _deltaX(_MAX_X / (DIM - 1))
  , _deltaY(_MAX_Y / (DIM - 1))
{
}

inline Visualization::~Visualization()
{
}

inline void Visualization::applyColorMap(const ColorMap* pColorMap, float value) const
{
  float R, G, B;
  pColorMap->map(value, R, G, B);
  glColor4f(R, G, B, pColorMap->getAlpha());
}

inline vec3f Visualization::crossproduct(const vec3f& vec1, const vec3f& vec2) const
{
  vec3f res;

  res[0] = (vec1[1] * vec2[2]) - (vec1[2] * vec2[1]);
  res[1] = (vec1[2] * vec2[0]) - (vec1[0] * vec2[2]);
  res[2] = (vec1[0] * vec2[1]) - (vec1[1] * vec2[0]);

  return res;
}

#endif /* VISUALIZATION_H_ */
