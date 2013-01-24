/*
 * GrDiffusionFTCS_Swap.h
 *
 *  Created on: Nov 23, 2010
 *      Author: pwolberg
 */

#ifndef GRDIFFUSIONFTCS_SWAP_H_
#define GRDIFFUSIONFTCS_SWAP_H_

#include "grdiffusion.h"

/**
* @brief Forward-Time Central-Space method for the heat equation [Reference](http://en.wikipedia.org/wiki/FTCS_scheme)
*/
class GrDiffusionFTCS_Swap: public GrDiffusion
{
private:
  /// Lower bound on when a chemical should be considered non-existent
  /// (to keep from going negative).
  const double _cutOffValue;

public:
  GrDiffusionFTCS_Swap();
  ~GrDiffusionFTCS_Swap();
  void diffuse(GrSimulationGrid& grid, const int time) const;
  DiffusionMethod getMethod() const;
  /*virtual*/
  GrDiffusion* clone() const
  {
    return new GrDiffusionFTCS_Swap(*this);
  }
};

inline DiffusionMethod GrDiffusionFTCS_Swap::getMethod() const
{
  return DIFF_REC_EQ_SWAP;
}

#endif /* GRDIFFUSIONFTCS_SWAP_H_ */
