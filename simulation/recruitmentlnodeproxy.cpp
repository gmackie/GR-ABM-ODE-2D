/*
 * recruitmentlnodeproxy.cpp
 *
 *  Created on: Mar 2, 2012
 *      Author: pwolberg
 */

#include "recruitmentlnodeproxy.h"
#include "stat.h"

using namespace std;

RecruitmentLnODEProxy::RecruitmentLnODEProxy()
: RecruitmentLnODE("", "")
{
	// TODO Auto-generated constructor stub
}

RecruitmentLnODEProxy::~RecruitmentLnODEProxy() {
	// TODO Auto-generated destructor stub
}

void RecruitmentLnODEProxy::solveODE(const int time, const Stats& statsPrevious, Stats& stats)
{
	//std::cerr << "_PARAM(PARAM_TCELL_LYMPH_PROXY_SCALING_START): " << _PARAM(PARAM_TCELL_LYMPH_PROXY_SCALING_START) << std::endl; //DBG
	//std::cerr << "_PARAM(PARAM_TCELL_LYMPH_PROXY_FACTOR_TIME_START): " << _PARAM(PARAM_TCELL_LYMPH_PROXY_FACTOR_TIME_START) << std::endl; //DBG
	//std::cerr << "_PARAM(PARAM_TCELL_LYMPH_PROXY_FACTOR_TIME_END): " << _PARAM(PARAM_TCELL_LYMPH_PROXY_FACTOR_TIME_END) << std::endl; //DBG
	//std::cerr << "_PARAM(PARAM_TCELL_LYMPH_PROXY_FACTOR_TIME_NONLINEAR_END): " << _PARAM(PARAM_TCELL_LYMPH_PROXY_FACTOR_TIME_NONLINEAR_END) << std::endl; //DBG
	//std::cerr << "_PARAM(PARAM_TCELL_LYMPH_PROXY_MTB_THRESHOLD): " << _PARAM(PARAM_TCELL_LYMPH_PROXY_MTB_THRESHOLD) << std::endl; //DBG

	Scalar totMtb = statsPrevious.getTotExtMtb() + statsPrevious.getTotIntMtb();
	//Scalar replicatingMtb =  totMtb - statsPrevious.getTotNonRepExtMtb();

	// In function RecruitmentLnODE::updateQueue the T gam flux is based on _odeInitialConditions[_idxEffectorT8] +  _odeInitialConditions[_idxEffectorTH1],
	// both of which are defined by the lymph node ODE. Here in the proxy calculation we calculate a single value for T gam flux,
	// so we store that in _odeInitialConditions[_idxEffectorT8] and set _odeInitialConditions[_idxEffectorTH1] to 0. Then in RecruitmentLnODE::updateQueue
	// for this proxy, the T gam flux is based on the single value calculated here.
	_odeInitialConditions[_idxEffectorTH1] = 0.0;
	if (time >= 0 && time < _PARAM(PARAM_TCELL_LYMPH_PROXY_SCALING_START))
	{
		_odeInitialConditions[_idxEffectorT8] = 0.0267 * totMtb;
		_odeInitialConditions[_idxCTL] = 0.0101 * totMtb;
		//std::cerr << "time: " << time << ", totMtb: " << totMtb << ", if1  _odeInitialConditions[_idxEffectorT8]: " << _odeInitialConditions[_idxEffectorT8] << std::endl; //DBG
	}
	else if (time >= _PARAM(PARAM_TCELL_LYMPH_PROXY_SCALING_START)
			&& time < _PARAM(PARAM_TCELL_LYMPH_PROXY_FACTOR_TIME_START))
	{
		if (totMtb >= _PARAM(PARAM_TCELL_LYMPH_PROXY_MTB_THRESHOLD))
		{
				_odeInitialConditions[_idxEffectorT8] = 0.04526 * totMtb + (0.0367 * time - 110);
				_odeInitialConditions[_idxCTL] = 0.017367 * totMtb + (0.0141 * time - 42.239);
				//std::cerr << "time: " << time << ", totMtb: " << totMtb << ", if2  _odeInitialConditions[_idxEffectorT8]: " << _odeInitialConditions[_idxEffectorT8] << std::endl; //DBG
		}
		else
		{
			_odeInitialConditions[_idxEffectorT8] = 0.0;
			_odeInitialConditions[_idxCTL] = 0.0;
			std::cout << "time: " << time << ", totMtb: " << totMtb << ", if5  _odeInitialConditions[_idxEffectorT8]: " << _odeInitialConditions[_idxEffectorT8] << std::endl; //DBG
		}
	}
	else if (time >= _PARAM(PARAM_TCELL_LYMPH_PROXY_FACTOR_TIME_START)
			&& time < _PARAM(PARAM_TCELL_LYMPH_PROXY_FACTOR_TIME_END))
	{
		if (totMtb >= _PARAM(PARAM_TCELL_LYMPH_PROXY_MTB_THRESHOLD))
		{
		Scalar timeDelta1 = _PARAM(PARAM_TCELL_LYMPH_PROXY_FACTOR_TIME_END) -
				_PARAM(PARAM_TCELL_LYMPH_PROXY_FACTOR_TIME_START);
		Scalar fluxFactor1 = (time - _PARAM(PARAM_TCELL_LYMPH_PROXY_FACTOR_TIME_START))/ timeDelta1;

		_odeInitialConditions[_idxEffectorT8] = fluxFactor1 * (0.337 * totMtb + 257.48);
		_odeInitialConditions[_idxCTL] = fluxFactor1 * (0.1302 * totMtb + 99.041);
		//std::cerr << "time: " << time << "fluxFactor1: " << fluxFactor1 << std::endl; //DBG
		//std::cerr << "time: " << time << ", totMtb: " << totMtb << ", if3  _odeInitialConditions[_idxEffectorT8]: " << _odeInitialConditions[_idxEffectorT8] << std::endl; //DBG
		}
		else
		{
			_odeInitialConditions[_idxEffectorT8] = 0.0;
			_odeInitialConditions[_idxCTL] = 0.0;
			std::cout << "time: " << time << ", totMtb: " << totMtb << ", if5  _odeInitialConditions[_idxEffectorT8]: " << _odeInitialConditions[_idxEffectorT8] << std::endl; //DBG
		}
	}
	else if (time >= _PARAM(PARAM_TCELL_LYMPH_PROXY_FACTOR_TIME_END)
			&& time < _PARAM(PARAM_TCELL_LYMPH_PROXY_FACTOR_TIME_NONLINEAR_END))
	{
		if (totMtb >= _PARAM(PARAM_TCELL_LYMPH_PROXY_MTB_THRESHOLD))
		{
			Scalar timeDelta2 = _PARAM(PARAM_TCELL_LYMPH_PROXY_FACTOR_TIME_NONLINEAR_END) -
					_PARAM(PARAM_TCELL_LYMPH_PROXY_FACTOR_TIME_END);
			Scalar fluxFactor2 = (time - _PARAM(PARAM_TCELL_LYMPH_PROXY_FACTOR_TIME_END))/ timeDelta2;
			//Scalar fluxFactor2 = std::min<Scalar>(0,(time - _PARAM(PARAM_TCELL_LYMPH_PROXY_FACTOR_TIME_END))/ timeDelta2);
			_odeInitialConditions[_idxEffectorT8] = fluxFactor2 * (0.0019 * pow(totMtb, 2) + 2.5275 * totMtb + 475.34);
			_odeInitialConditions[_idxCTL] = fluxFactor2 * (0.0007 * pow(totMtb, 2) + 0.9734 * totMtb + 182.98);
			//std::cerr << "time: " << time << "fluxFactor2: " << fluxFactor2 << std::endl; //DBG
			std::cerr << "time: " << time << ", totMtb: " << totMtb << ", if4  _odeInitialConditions[_idxEffectorT8]: " << _odeInitialConditions[_idxEffectorT8] << std::endl; //DBG
		}
		else
		{
			_odeInitialConditions[_idxEffectorT8] = 0.0;
			_odeInitialConditions[_idxCTL] = 0.0;
			std::cout << "time: " << time << ", totMtb: " << totMtb << ", if5  _odeInitialConditions[_idxEffectorT8]: " << _odeInitialConditions[_idxEffectorT8] << std::endl; //DBG
		}
	}
	else
	{
		if (totMtb >= _PARAM(PARAM_TCELL_LYMPH_PROXY_MTB_THRESHOLD))
		{
			_odeInitialConditions[_idxEffectorT8] = (0.0019 * pow(totMtb, 2) + 2.5275 * totMtb + 475.34);
			_odeInitialConditions[_idxCTL] = (0.0007 * pow(totMtb, 2) + 0.9734 * totMtb + 182.98);
			std::cout << "time: " << time << ", totMtb: " << totMtb << ", if6  _odeInitialConditions[_idxEffectorT8]: " << _odeInitialConditions[_idxEffectorT8] << std::endl; //DBG
		}
		else
		{
			_odeInitialConditions[_idxEffectorT8] = 0.0;
			_odeInitialConditions[_idxCTL] = 0.0;
			std::cout << "time: " << time << ", totMtb: " << totMtb << ", if7  _odeInitialConditions[_idxEffectorT8]: " << _odeInitialConditions[_idxEffectorT8] << std::endl; //DBG
		}
	}
    stats.setTH1lung(_odeInitialConditions[_idxEffectorTH1]);
    stats.setT8lung(_odeInitialConditions[_idxEffectorT8]);
    stats.setTClung(_odeInitialConditions[_idxCTL]);
}
