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

class GrDiffusion
{
public:
	GrDiffusion();
	virtual ~GrDiffusion();
	virtual void diffuse(GrSimulationGrid& grid) const = 0;
	virtual DiffusionMethod getMethod() const = 0;
};

#endif /* GRDIFFUSION_H_ */
