/*
 * recruitmentprob.h
 *
 *  Created on: Jun 21, 2010
 *      Author: mohammed
 */

#ifndef RECRUITMENTPROB_H_
#define RECRUITMENTPROB_H_

#include "gr.h"
#include "grgrid.h"
#include "grstat.h"
#include "recruitmentbase.h"

class RecruitmentProb : public RecruitmentBase
{
private:
	void recruitMac(GrSimulation& sim, GridCell* pSource);
	void recruitTcell(GrSimulation& sim, GridCell* pSource);

public:
	RecruitmentProb();
	virtual ~RecruitmentProb();
	void recruit(GrSimulation& sim);
};

#endif /* RECRUITMENTPROB_H_ */
