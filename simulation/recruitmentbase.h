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

	static bool MacRecruitmentThreshold(const GrGrid& g, const Pos& pSource);
	static bool TgamRecruitmentThreshold(const GrGrid& g, const Pos& pSource);
	static bool TcytRecruitmentThreshold(const GrGrid& g, const Pos& pSource);
	static bool TregRecruitmentThreshold(const GrGrid& g, const Pos& pSource);

	static bool MacRecruitmentThreshold(const GrGrid& g, const Pos& pSource, double& rThreshold);
	static bool TgamRecruitmentThreshold(const GrGrid& g, const Pos& pSource, double& rThreshold);
	static bool TcytRecruitmentThreshold(const GrGrid& g, const Pos& pSource, double& rThreshold);
	static bool TregRecruitmentThreshold(const GrGrid& g, const Pos& pSource, double& rThreshold);
};

inline bool RecruitmentBase::MacRecruitmentThreshold(const GrGrid& grid, const Pos& pSource, double& rThreshold)
{
	rThreshold = _PARAM(PARAM_GR_WEIGHT_TNF_RECRUITMENT) * grid.TNF(pSource) +
			_PARAM(PARAM_GR_WEIGHT_CCL2_RECRUITMENT) * grid.CCL2(pSource) +
			_PARAM(PARAM_GR_WEIGHT_CCL5_RECRUITMENT) * grid.CCL5(pSource);

	return rThreshold >= _PARAM(PARAM_MAC_THRESHOLD_RECRUITMENT);
}

inline bool RecruitmentBase::TgamRecruitmentThreshold(const GrGrid& grid, const Pos& pSource, double& rThreshold)
{
	rThreshold = _PARAM(PARAM_GR_WEIGHT_TNF_RECRUITMENT) * grid.TNF(pSource) +
			_PARAM(PARAM_GR_WEIGHT_CCL2_RECRUITMENT) * grid.CCL2(pSource) +
			_PARAM(PARAM_GR_WEIGHT_CCL5_RECRUITMENT) * grid.CCL5(pSource) +
			_PARAM(PARAM_GR_WEIGHT_CXCL9_RECRUITMENT) * grid.CXCL9(pSource);

	return rThreshold >= _PARAM(PARAM_TGAM_THRESHOLD_RECRUITMENT);
}

inline bool RecruitmentBase::TcytRecruitmentThreshold(const GrGrid& grid, const Pos& pSource, double& rThreshold)
{
	rThreshold = _PARAM(PARAM_GR_WEIGHT_TNF_RECRUITMENT) * grid.TNF(pSource) +
		//_PARAM(PARAM_GR_WEIGHT_CCL2_RECRUITMENT) * grid.CCL2() +
		_PARAM(PARAM_GR_WEIGHT_CCL5_RECRUITMENT) * grid.CCL5(pSource) +
		_PARAM(PARAM_GR_WEIGHT_CXCL9_RECRUITMENT) * grid.CXCL9(pSource);

	return rThreshold >= _PARAM(PARAM_TCYT_THRESHOLD_RECRUITMENT);
}

inline bool RecruitmentBase::TregRecruitmentThreshold(const GrGrid& grid, const Pos& pSource, double& rThreshold)
{
	rThreshold = _PARAM(PARAM_GR_WEIGHT_TNF_RECRUITMENT) * grid.TNF(pSource) +
		_PARAM(PARAM_GR_WEIGHT_CCL5_RECRUITMENT) * grid.CCL5(pSource);

	return rThreshold >= _PARAM(PARAM_TREG_THRESHOLD_RECRUITMENT);
}

inline bool RecruitmentBase::MacRecruitmentThreshold(const GrGrid& grid, const Pos& pSource)
{
	double threshold;
	return MacRecruitmentThreshold(grid, pSource, threshold);
}

inline bool RecruitmentBase::TgamRecruitmentThreshold(const GrGrid& grid, const Pos& pSource)
{
	double threshold;
	return TgamRecruitmentThreshold(grid, pSource, threshold);
}

inline bool RecruitmentBase::TcytRecruitmentThreshold(const GrGrid& grid, const Pos& pSource)
{
	double threshold;
	return TcytRecruitmentThreshold(grid, pSource, threshold);
}

inline bool RecruitmentBase::TregRecruitmentThreshold(const GrGrid& grid, const Pos& pSource)
{
	double threshold;
	return TregRecruitmentThreshold(grid, pSource, threshold);
}

#endif /* RECRUITMENTBASE_H_ */
