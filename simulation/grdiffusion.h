/*
 * grdiffusion.h
 *
 *  Created on: Feb 2, 2010
 *      Author: mohammed
 */

#ifndef GRDIFFUSION_H_
#define GRDIFFUSION_H_

#include "gr.h"
#include "grgrid.h"

class GrDiffusion
{
public:
	GrDiffusion();
	virtual ~GrDiffusion();
	virtual void diffuse(GrGrid& grid) const = 0;
	virtual DiffusionMethod getMethod() const = 0;
};

#endif /* GRDIFFUSION_H_ */
