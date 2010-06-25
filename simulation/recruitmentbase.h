/*
 * recruitmentbase.h
 *
 *  Created on: Jun 21, 2010
 *      Author: mohammed
 */

#ifndef RECRUITMENTBASE_H_
#define RECRUITMENTBASE_H_

#include "gr.h"
#include "grgrid.h"
#include "grstat.h"
#include "grsimulation.h"

class RecruitmentBase
{
public:
	RecruitmentBase();
	virtual ~RecruitmentBase();
	virtual void recruit(GrSimulation& sim) = 0;

	static bool MacRecruitmentThreshold(const GridCell* pSource);
	static bool TgamRecruitmentThreshold(const GridCell* pSource);
	static bool TcytRecruitmentThreshold(const GridCell* pSource);
	static bool TregRecruitmentThreshold(const GridCell* pSource);

	static bool MacRecruitmentThreshold(const GridCell* pSource, double& rThreshold);
	static bool TgamRecruitmentThreshold(const GridCell* pSource, double& rThreshold);
	static bool TcytRecruitmentThreshold(const GridCell* pSource, double& rThreshold);
	static bool TregRecruitmentThreshold(const GridCell* pSource, double& rThreshold);
};

inline bool RecruitmentBase::MacRecruitmentThreshold(const GridCell* pSource, double& rThreshold)
{
	rThreshold = _PARAM(PARAM_GR_WEIGHT_TNF_RECRUITMENT) * pSource->getTNF() +
			_PARAM(PARAM_GR_WEIGHT_CCL2_RECRUITMENT) * pSource->getCCL2() +
			_PARAM(PARAM_GR_WEIGHT_CCL5_RECRUITMENT) * pSource->getCCL5();

	return rThreshold >= _PARAM(PARAM_MAC_THRESHOLD_RECRUITMENT);
}

inline bool RecruitmentBase::TgamRecruitmentThreshold(const GridCell* pSource, double& rThreshold)
{
	rThreshold = _PARAM(PARAM_GR_WEIGHT_TNF_RECRUITMENT) * pSource->getTNF() +
			_PARAM(PARAM_GR_WEIGHT_CCL2_RECRUITMENT) * pSource->getCCL2() +
			_PARAM(PARAM_GR_WEIGHT_CCL5_RECRUITMENT) * pSource->getCCL5() +
			_PARAM(PARAM_GR_WEIGHT_CXCL9_RECRUITMENT) * pSource->getCXCL9();

	return rThreshold >= _PARAM(PARAM_TGAM_THRESHOLD_RECRUITMENT);
}

inline bool RecruitmentBase::TcytRecruitmentThreshold(const GridCell* pSource, double& rThreshold)
{
	rThreshold = _PARAM(PARAM_GR_WEIGHT_TNF_RECRUITMENT) * pSource->getTNF() +
		//_PARAM(PARAM_GR_WEIGHT_CCL2_RECRUITMENT) * pSource->getCCL2() +
		_PARAM(PARAM_GR_WEIGHT_CCL5_RECRUITMENT) * pSource->getCCL5() +
		_PARAM(PARAM_GR_WEIGHT_CXCL9_RECRUITMENT) * pSource->getCXCL9();

	return rThreshold >= _PARAM(PARAM_TCYT_THRESHOLD_RECRUITMENT);
}

inline bool RecruitmentBase::TregRecruitmentThreshold(const GridCell* pSource, double& rThreshold)
{
	rThreshold = _PARAM(PARAM_GR_WEIGHT_TNF_RECRUITMENT) * pSource->getTNF() +
		_PARAM(PARAM_GR_WEIGHT_CCL5_RECRUITMENT) * pSource->getCCL5();

	return rThreshold >= _PARAM(PARAM_TREG_THRESHOLD_RECRUITMENT);
}

inline bool RecruitmentBase::MacRecruitmentThreshold(const GridCell* pSource)
{
	double threshold;
	return MacRecruitmentThreshold(pSource, threshold);
}

inline bool RecruitmentBase::TgamRecruitmentThreshold(const GridCell* pSource)
{
	double threshold;
	return TgamRecruitmentThreshold(pSource, threshold);
}

inline bool RecruitmentBase::TcytRecruitmentThreshold(const GridCell* pSource)
{
	double threshold;
	return TcytRecruitmentThreshold(pSource, threshold);
}

inline bool RecruitmentBase::TregRecruitmentThreshold(const GridCell* pSource)
{
	double threshold;
	return TregRecruitmentThreshold(pSource, threshold);
}

#endif /* RECRUITMENTBASE_H_ */
