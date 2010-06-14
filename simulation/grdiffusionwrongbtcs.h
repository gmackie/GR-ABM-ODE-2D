/*
 * grdiffusionwrongbtcs.h
 *
 *  Created on: Feb 2, 2010
 *      Author: mohammed
 */

#ifndef GRDIFFUSIONWRONGBTCS_H_
#define GRDIFFUSIONWRONGBTCS_H_

#include "grdiffusion.h"

class GrDiffusionWrongBTCS: public GrDiffusion
{
public:
	GrDiffusionWrongBTCS();
	virtual ~GrDiffusionWrongBTCS();
	void diffuse(GrGrid& grid) const;
	DiffusionMethod getMethod() const;
};

inline DiffusionMethod GrDiffusionWrongBTCS::getMethod() const
{
	return DIFF_SOR_WRONG;
}

#endif /* GRDIFFUSIONWRONGBTCS_H_ */
