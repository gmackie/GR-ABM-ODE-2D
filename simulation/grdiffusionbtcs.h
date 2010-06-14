/*
 * grdiffusionbtcs.h
 *
 *  Created on: Feb 2, 2010
 *      Author: mohammed
 */

#ifndef GRDIFFUSIONBTCS_H_
#define GRDIFFUSIONBTCS_H_

#include "grdiffusion.h"

class GrDiffusionBTCS: public GrDiffusion
{
public:
	GrDiffusionBTCS();
	virtual ~GrDiffusionBTCS();
	void diffuse(GrGrid& grid) const;
	DiffusionMethod getMethod() const;
};

inline DiffusionMethod GrDiffusionBTCS::getMethod() const
{
	return DIFF_SOR_CORRECT;
}

#endif /* GRDIFFUSIONBTCS_H_ */
