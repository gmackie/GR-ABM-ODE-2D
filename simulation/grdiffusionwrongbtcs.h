/*
 * grdiffusionwrongbtcs.h
 *
 *  Created on: Feb 2, 2010
 *      Author: mohammed
 */

#ifndef GRDIFFUSIONWRONGBTCS_H_
#define GRDIFFUSIONWRONGBTCS_H_

#include "grdiffusion.h"

/**
* @brief Do not use this method, it is not updated, nor was it
* mathematically correct to begin with
* @deprecated
*/
class GrDiffusionWrongBTCS: public GrDiffusion
{
public:
  GrDiffusionWrongBTCS();
  virtual ~GrDiffusionWrongBTCS();
  void diffuse(GrSimulationGrid& grid, const int time) const;
  DiffusionMethod getMethod() const;
  GrDiffusion* clone() const
  {
    return new GrDiffusionWrongBTCS(*this);
  }
};

inline DiffusionMethod GrDiffusionWrongBTCS::getMethod() const
{
  return DIFF_SOR_WRONG;
}

#endif /* GRDIFFUSIONWRONGBTCS_H_ */
