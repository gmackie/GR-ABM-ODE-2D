/*
 * GrDiffusionFTCS_Swap.h
 *
 *  Created on: Nov 23, 2010
 *      Author: pwolberg
 */

#ifndef GRDIFFUSIONFTCS_SWAP_H_
#define GRDIFFUSIONFTCS_SWAP_H_

#include "grdiffusion.h"

class GrDiffusionFTCS_Swap: public GrDiffusion
{
private:
	const double _cutOffValue;

public:
	GrDiffusionFTCS_Swap();
	virtual ~GrDiffusionFTCS_Swap();
	void diffuse(GrSimulationGrid& grid) const;
	DiffusionMethod getMethod() const;
};

inline DiffusionMethod GrDiffusionFTCS_Swap::getMethod() const
{
	return DIFF_REC_EQ_SWAP;
}

#endif /* GRDIFFUSIONFTCS_SWAP_H_ */
