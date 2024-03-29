/*
 * grdiffusion.h
 *
 *  Created on: Feb 2, 2010
 *      Author: mohammed
 */

#ifndef GRDIFFUSION_H_
#define GRDIFFUSION_H_

#include "gr.h"
#include "grsimulationgrid.h"

/**
* @brief Diffusion base class for defining diffusion methods
* @todo Redefine this class to deal with individual chemicals, abstracting from
* the diffusing chemicals (let GrSimulation or GrGrid handle it).  Would be best
* to have the diffusion method take in two ScalarGrid, one read one write, and
* work from there.
*/
class GrDiffusion
{
public:
  GrDiffusion();
  virtual ~GrDiffusion();
  /**
  * @brief Diffuse selected chemicals on the grid
  */
  virtual void diffuse(GrSimulationGrid& grid, const int time) const = 0;
  /**
  * @return The enumerated type of the method
  */
  virtual DiffusionMethod getMethod() const = 0;
  /**
  * @brief Deep copy of the diffusion method
  */
  virtual GrDiffusion* clone() const = 0;
};

#endif /* GRDIFFUSION_H_ */
