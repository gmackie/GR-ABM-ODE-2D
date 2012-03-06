/*
 * recruitmentlnodeproxy.cpp
 *
 *  Created on: Mar 2, 2012
 *      Author: pwolberg
 */

#include "recruitmentlnodeproxy.h"

RecruitmentLnODEProxy::RecruitmentLnODEProxy()
: RecruitmentLnODE("", "")
{
	// TODO Auto-generated constructor stub
}

RecruitmentLnODEProxy::~RecruitmentLnODEProxy() {
	// TODO Auto-generated destructor stub
}

void RecruitmentLnODEProxy::solveODE(const int time, GrStat& stats)
{
	int boundaryTime = 20 * TIME_STEPS_PER_DAY;

	FLOAT_TYPE totMtb = stats. getTotExtMtb() + stats.getTotIntMtb();
	FLOAT_TYPE replicatingMtb =  totMtb - stats.getTotNonRepExtMtb();

	// In function RecruitmentLnODE::updateQueue the T gam flux is based on _odeInitialConditions[_idxEffectorT8] +  _odeInitialConditions[_idxEffectorTH1],
	// both of which are defined by the lymph node ODE. Here in the proxy calculation we calculate a single value for T gam flux,
	// so we store that in _odeInitialConditions[_idxEffectorT8] and set _odeInitialConditions[_idxEffectorTH1] to 0. Then in RecruitmentLnODE::updateQueue
	// for this proxy, the T gam flux is based in the single value calculated here.
	_odeInitialConditions[_idxEffectorTH1] = 0.0;
	if (time >= 0 && time < boundaryTime)
	{
		if (replicatingMtb > 0.0)
		{
			_odeInitialConditions[_idxEffectorT8] = std::max<FLOAT_TYPE>(0, (2.0496 * log(replicatingMtb) - 6.5451));
			_odeInitialConditions[_idxCTL] = std::max<FLOAT_TYPE>(0, (0.779 * log(replicatingMtb) - 2.491));
		}
		else
		{
			_odeInitialConditions[_idxEffectorT8] = 0.0;
			_odeInitialConditions[_idxCTL] = 0.0;
		}
	}
	else
	{
		_odeInitialConditions[_idxEffectorT8] = 0.0012 * pow(replicatingMtb, 2) + 2.0676 * replicatingMtb + 144.54;
		_odeInitialConditions[_idxCTL] = 0.0005 * pow(replicatingMtb, 2) + 0.7956 * replicatingMtb + 55.803;
	}
}
