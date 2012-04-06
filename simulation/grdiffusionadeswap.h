/*
 * GrDiffusionADE_Swap.h
 *
 *  Created on: Mar 8, 2012
 *      Author: ncilfone
 */

#ifndef GRDIFFUSIONADE_SWAP_H_
#define GRDIFFUSIONADE_SWAP_H_

#include "grdiffusion.h"

class GrDiffusionADE_Swap: public GrDiffusion
{
private:
	const double _cutOffValue;

public:
	GrDiffusionADE_Swap();
	virtual ~GrDiffusionADE_Swap();
	void diffuse(GrSimulationGrid& grid) const;
	DiffusionMethod getMethod() const;
};

inline DiffusionMethod GrDiffusionADE_Swap::getMethod() const
{
	return DIFF_ADE_SWAP;
}


#endif /* GRDIFFUSIONADE_SWAP_H_ */
