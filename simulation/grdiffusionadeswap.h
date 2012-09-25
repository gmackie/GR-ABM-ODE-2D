/*
 * GrDiffusionADE_Swap.h
 *
 *  Created on: Mar 8, 2012
 *      Author: ncilfone
 */

#ifndef GRDIFFUSIONADE_SWAP_H_
#define GRDIFFUSIONADE_SWAP_H_

#include "grdiffusion.h"

/**
* @brief 2-D Alternating Direction Explicit Algorithm [Reference](about:blank)
* @details For reference on the alogrithm, read the following paper:
* On the Solution of the Diffusion Equations by Numerical Methods
* H.Z. BARAKAT and J.A. Clark
* Journal of Heat Transfer November, 1966 pp 421-427
* @note
* Note that the derivation in this paper uses i=cols and j=rows
* We use i=rows and j=cols thus our coefficients b and c are different
* than what is shown in the paper (they are switched)
*/
class GrDiffusionADE_Swap: public GrDiffusion
{
private:
  /// Lower bound on when a chemical should be considered non-existent
  /// (to keep from going negative).
  const double _cutOffValue;

public:
  GrDiffusionADE_Swap();
  virtual ~GrDiffusionADE_Swap();
  void diffuse(GrSimulationGrid& grid) const;
  /*virtual*/
  GrDiffusion* clone() const
  {
    return new GrDiffusionADE_Swap(*this);
  }
  DiffusionMethod getMethod() const;
};

inline DiffusionMethod GrDiffusionADE_Swap::getMethod() const
{
  return DIFF_ADE_SWAP;
}


#endif /* GRDIFFUSIONADE_SWAP_H_ */
