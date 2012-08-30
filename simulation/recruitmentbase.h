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

using namespace std;

class RecruitmentBase
{
public:
	RecruitmentBase();
	virtual ~RecruitmentBase();

	virtual RecruitmentMethod getMethod() const = 0;
	virtual void serialize(std::ostream& out) const = 0;
	virtual void deserialize(std::istream& in) = 0;

    virtual void recruit(GrSimulation& sim, int time) = 0;
  virtual RecruitmentBase* clone() const = 0;
    
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
    
    double intpart1, intpart2, Store1, Store2, StorePower;
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
        
        modf(Store1, &intpart1);
        modf(Store2, &intpart2);
        
        intpart1Store = (int)intpart1;
        intpart2Store = (int)intpart2;
        
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
//    double Vmax = (1.0)/(3.0); // Temporary Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1

    rThreshold = _PARAM(PARAM_GR_WEIGHT_TNF_RECRUITMENT) * grid.TNF(pSource) +
            _PARAM(PARAM_GR_WEIGHT_CCL2_RECRUITMENT) * grid.CCL2(pSource) +
            _PARAM(PARAM_GR_WEIGHT_CCL5_RECRUITMENT) * grid.CCL5(pSource);
    
//    std::cout << "Mac Recruitment Threshold: " << rThreshold << std::endl;

//    rThreshold = Vmax * ((grid.CCL2(pSource)/(grid.CCL2(pSource) + _PARAM(PARAM_MAC_HALF_SAT_CHEMOKINE))) +
//                         (grid.CCL5(pSource)/(grid.CCL5(pSource) + _PARAM(PARAM_MAC_HALF_SAT_CHEMOKINE))) +
//                         (grid.TNF(pSource)/(grid.TNF(pSource) + _PARAM(PARAM_MAC_HALF_SAT_TNF))));


	return intCompareGTEQ(rThreshold, _PARAM(PARAM_MAC_THRESHOLD_RECRUITMENT));
}

inline bool RecruitmentBase::TgamRecruitmentThreshold(const GrGrid& grid, const Pos& pSource, double& rThreshold)
{
//    double Vmax = (1.0)/(4.0); // Temporary Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1

    rThreshold = _PARAM(PARAM_GR_WEIGHT_TNF_RECRUITMENT) * grid.TNF(pSource) +
            _PARAM(PARAM_GR_WEIGHT_CCL2_RECRUITMENT) * grid.CCL2(pSource) +
            _PARAM(PARAM_GR_WEIGHT_CCL5_RECRUITMENT) * grid.CCL5(pSource) +
            _PARAM(PARAM_GR_WEIGHT_CXCL9_RECRUITMENT) * grid.CXCL9(pSource);

//     std::cout << "Tgam Recruitment Threshold: " << rThreshold << std::endl;

/*    rThreshold = Vmax * ((grid.CCL2(pSource)/(grid.CCL2(pSource) + _PARAM(PARAM_MAC_HALF_SAT_CHEMOKINE))) +
                         (grid.CCL5(pSource)/(grid.CCL5(pSource) + _PARAM(PARAM_MAC_HALF_SAT_CHEMOKINE))) +
                         (grid.TNF(pSource)/(grid.TNF(pSource) + _PARAM(PARAM_MAC_HALF_SAT_TNF))) +
                          (grid.CXCL9(pSource)/(grid.CXCL9(pSource) + _PARAM(PARAM_MAC_HALF_SAT_CHEMOKINE))))*/;


	return intCompareGTEQ(rThreshold, _PARAM(PARAM_TGAM_THRESHOLD_RECRUITMENT));
}

inline bool RecruitmentBase::TcytRecruitmentThreshold(const GrGrid& grid, const Pos& pSource, double& rThreshold)
{
//    double Vmax = (1.0)/(3.0); // Temporary Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1

    rThreshold = _PARAM(PARAM_GR_WEIGHT_TNF_RECRUITMENT) * grid.TNF(pSource) +
        //_PARAM(PARAM_GR_WEIGHT_CCL2_RECRUITMENT) * grid.CCL2() +
        _PARAM(PARAM_GR_WEIGHT_CCL5_RECRUITMENT) * grid.CCL5(pSource) +
        _PARAM(PARAM_GR_WEIGHT_CXCL9_RECRUITMENT) * grid.CXCL9(pSource);


//    rThreshold = Vmax * ((grid.CCL5(pSource)/(grid.CCL5(pSource) + _PARAM(PARAM_MAC_HALF_SAT_CHEMOKINE))) +
//                         (grid.TNF(pSource)/(grid.TNF(pSource) + _PARAM(PARAM_MAC_HALF_SAT_TNF))) +
//                          (grid.CXCL9(pSource)/(grid.CXCL9(pSource) + _PARAM(PARAM_MAC_HALF_SAT_CHEMOKINE))));

	return intCompareGTEQ(rThreshold, _PARAM(PARAM_TCYT_THRESHOLD_RECRUITMENT));
}

inline bool RecruitmentBase::TregRecruitmentThreshold(const GrGrid& grid, const Pos& pSource, double& rThreshold)
{
//    double Vmax = (1.0)/(2.0); // Temporary Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1

    rThreshold = _PARAM(PARAM_GR_WEIGHT_TNF_RECRUITMENT) * grid.TNF(pSource) +
        _PARAM(PARAM_GR_WEIGHT_CCL5_RECRUITMENT) * grid.CCL5(pSource);


//    rThreshold = Vmax * ((grid.CCL5(pSource)/(grid.CCL5(pSource) + _PARAM(PARAM_MAC_HALF_SAT_CHEMOKINE))) +
//                         (grid.TNF(pSource)/(grid.TNF(pSource) + _PARAM(PARAM_MAC_HALF_SAT_TNF))));

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
