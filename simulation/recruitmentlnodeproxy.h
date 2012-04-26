/*
 * recruitmentlnodeproxy.h
 *
 *  Created on: Mar 2, 2012
 *      Author: pwolberg
 */

#include "recruitmentlnode.h"

#ifndef RECRUITMENTLNODEPROXY_H_
#define RECRUITMENTLNODEPROXY_H_

class RecruitmentLnODEProxy  : public RecruitmentLnODE
{
private:
	virtual void solveODE(const int time, const Stats& statsPrevious, Stats& stats);

public:
	RecruitmentLnODEProxy();
	virtual ~RecruitmentLnODEProxy();
};

#endif /* RECRUITMENTLNODEPROXY_H_ */
