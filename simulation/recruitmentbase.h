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
#include "stat.h"
#include "grsimulation.h"

using namespace std;

class RecruitmentBase
{
public:
	RecruitmentBase();
	virtual ~RecruitmentBase();
	virtual void recruit(GrSimulation& sim) = 0;
    
    static bool intCompareGTEQ(const double param1, const double param2);

	static bool MacRecruitmentThreshold(const GrGrid& g, const Pos& pSource);
	static bool TgamRecruitmentThreshold(const GrGrid& g, const Pos& pSource);
	static bool TcytRecruitmentThreshold(const GrGrid& g, const Pos& pSource);
	static bool TregRecruitmentThreshold(const GrGrid& g, const Pos& pSource);

	static bool MacRecruitmentThreshold(const GrGrid& g, const Pos& pSource, double& rThreshold);
	static bool TgamRecruitmentThreshold(const GrGrid& g, const Pos& pSource, double& rThreshold);
	static bool TcytRecruitmentThreshold(const GrGrid& g, const Pos& pSource, double& rThreshold);
	static bool TregRecruitmentThreshold(const GrGrid& g, const Pos& pSource, double& rThreshold);
};


inline bool RecruitmentBase::intCompareGTEQ(const double param1, const double param2)
{
    // Compares two values based on the number of sig figs we hold in gr.h (ABS_TOL)
    // If the values are not within 2 orders of magnitude we do not convert to ints
    // since it should not matter
    // COMPARES PARAM1 >= PARAM2 and returns bool based on this evaluation
    
    double intpart1, intpart2, fracpart1, fracpart2, Store1, Store2, StorePower;
    int intpart1Store, intpart2Store;
    
    bool result = 0;
    
    if (fabs(floor(log10(param1)) - floor(log10(param2))) < 2  ) 
    {
        if (floor(log10(param1)) < floor(log10(param2)))
        {
            StorePower = floor(log10(param1));
        }
        else
        {
            StorePower = floor(log10(param2));
        }
        
        Store1 = (param1 * (ABS_TOL/(pow(10,StorePower))));
        Store2 = (param2 * (ABS_TOL/(pow(10,StorePower))));
        
        fracpart1 = modf(Store1, &intpart1);
        fracpart2 = modf(Store2, &intpart2);
        
        intpart1Store = intpart1;
        intpart2Store = intpart2;
        
        if (intpart1Store >= intpart2Store) 
        {
            result = 1;
        }
    }
    else
    {
        if (param1 >= param2) 
        {
            result = 1;
        }
    }
    return result;
}




inline bool RecruitmentBase::MacRecruitmentThreshold(const GrGrid& grid, const Pos& pSource, double& rThreshold)
{
	rThreshold = _PARAM(PARAM_GR_WEIGHT_TNF_RECRUITMENT) * grid.TNF(pSource) +
			_PARAM(PARAM_GR_WEIGHT_CCL2_RECRUITMENT) * grid.CCL2(pSource) +
			_PARAM(PARAM_GR_WEIGHT_CCL5_RECRUITMENT) * grid.CCL5(pSource);
    
	return intCompareGTEQ(rThreshold, _PARAM(PARAM_MAC_THRESHOLD_RECRUITMENT));
}

inline bool RecruitmentBase::TgamRecruitmentThreshold(const GrGrid& grid, const Pos& pSource, double& rThreshold)
{
	rThreshold = _PARAM(PARAM_GR_WEIGHT_TNF_RECRUITMENT) * grid.TNF(pSource) +
			_PARAM(PARAM_GR_WEIGHT_CCL2_RECRUITMENT) * grid.CCL2(pSource) +
			_PARAM(PARAM_GR_WEIGHT_CCL5_RECRUITMENT) * grid.CCL5(pSource) +
			_PARAM(PARAM_GR_WEIGHT_CXCL9_RECRUITMENT) * grid.CXCL9(pSource);

	return intCompareGTEQ(rThreshold, _PARAM(PARAM_TGAM_THRESHOLD_RECRUITMENT));
}

inline bool RecruitmentBase::TcytRecruitmentThreshold(const GrGrid& grid, const Pos& pSource, double& rThreshold)
{
	rThreshold = _PARAM(PARAM_GR_WEIGHT_TNF_RECRUITMENT) * grid.TNF(pSource) +
		//_PARAM(PARAM_GR_WEIGHT_CCL2_RECRUITMENT) * grid.CCL2() +
		_PARAM(PARAM_GR_WEIGHT_CCL5_RECRUITMENT) * grid.CCL5(pSource) +
		_PARAM(PARAM_GR_WEIGHT_CXCL9_RECRUITMENT) * grid.CXCL9(pSource);

	return intCompareGTEQ(rThreshold, _PARAM(PARAM_TCYT_THRESHOLD_RECRUITMENT));
}

inline bool RecruitmentBase::TregRecruitmentThreshold(const GrGrid& grid, const Pos& pSource, double& rThreshold)
{
	rThreshold = _PARAM(PARAM_GR_WEIGHT_TNF_RECRUITMENT) * grid.TNF(pSource) +
		_PARAM(PARAM_GR_WEIGHT_CCL5_RECRUITMENT) * grid.CCL5(pSource);

	return intCompareGTEQ(rThreshold, _PARAM(PARAM_TREG_THRESHOLD_RECRUITMENT));
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
