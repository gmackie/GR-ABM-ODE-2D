/*
 * grdiffusionftcs.h
 *
 *  Created on: Feb 2, 2010
 *      Author: mohammed
 */

#ifndef GRDIFFUSIONFTCS_H_
#define GRDIFFUSIONFTCS_H_

#include "grdiffusion.h"
#include "grgrid.h"

/**
 * Computes diffusion using Forward Time Centered Space
 *
 * _mu = D * dt / (dx)^2, where dt is 6 s and dx = 20 \mu m and D is in cm^2/s
 */
class GrDiffusionFTCS: public GrDiffusion
{
private:
	const double _cutOffValue;
	mutable GrGrid _cpyGrid;

public:
	GrDiffusionFTCS();
	virtual ~GrDiffusionFTCS();
	void diffuse(GrGrid& grid) const;
	DiffusionMethod getMethod() const;
};

inline DiffusionMethod GrDiffusionFTCS::getMethod() const
{
	return DIFF_REC_EQ;
}

#endif /* GRDIFFUSIONFTCS_H_ */
