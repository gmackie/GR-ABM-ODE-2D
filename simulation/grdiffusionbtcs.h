/*
 * grdiffusionbtcs.h
 *
 *  Created on: Feb 2, 2010
 *      Author: mohammed
 */

#ifndef GRDIFFUSIONBTCS_H_
#define GRDIFFUSIONBTCS_H_

#include "grdiffusion.h"

/**
* @brief Do not use this method, it is not updated, nor was it
* mathematically correct to begin with
* @deprecated
*/
class GrDiffusionBTCS: public GrDiffusion
{
public:
  GrDiffusionBTCS();
  virtual ~GrDiffusionBTCS();
  void diffuse(GrSimulationGrid& grid, const int time) const;
  DiffusionMethod getMethod() const;
  GrDiffusion* clone() const
  {
    return new GrDiffusionBTCS(*this);
  }
};

inline DiffusionMethod GrDiffusionBTCS::getMethod() const
{
  return DIFF_SOR_CORRECT;
}

#endif /* GRDIFFUSIONBTCS_H_ */
