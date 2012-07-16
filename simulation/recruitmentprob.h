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
#include "stat.h"
#include "recruitmentbase.h"

class RecruitmentProb : public RecruitmentBase
{
private:
	void recruitMac(GrSimulation& sim, const Pos& pSource);
	void recruitTcell(GrSimulation& sim, const Pos& pSource);
    void recruitCellMac(GrSimulation& sim, const Pos& pSource);
    void recruitCellTgam(GrSimulation& sim, const Pos& pSource);
    void recruitCellTcyt(GrSimulation& sim, const Pos& pSource);
    void recruitCellTreg(GrSimulation& sim, const Pos& pSource);

public:
	RecruitmentProb();
	virtual ~RecruitmentProb();
    void recruit(GrSimulation& sim);
    void recruit(GrSimulation& sim, const Pos& pSource);

    static bool intCompareGTEQ(const double param1, const double param2);
    static bool MacThresholdRecNew(const GrGrid& grid, const Pos& pSource, double& rThreshold);
    static bool TgamThresholdRecNew(const GrGrid& grid, const Pos& pSource, double& rThreshold);
    static bool TcytThresholdRecNew(const GrGrid& grid, const Pos& pSource, double& rThreshold);
    static bool TregThresholdRecNew(const GrGrid& grid, const Pos& pSource, double& rThreshold);
	RecruitmentMethod getMethod() const;
	void serialize(std::ostream&) const;
	void deserialize(std::istream&);

    static bool MacThresholdRecNew(const GrGrid& grid, const Pos& pSource);
    static bool TgamThresholdRecNew(const GrGrid& grid, const Pos& pSource);
    static bool TcytThresholdRecNew(const GrGrid& grid, const Pos& pSource);
    static bool TregThresholdRecNew(const GrGrid& grid, const Pos& pSource);

  /*virtual*/ RecruitmentBase* clone() const { return new RecruitmentProb(*this); }
protected:

    bool PossibleRecruitMac(GrSimulation& sim, const Pos& pSource);
    bool PossibleRecruitTgam(GrSimulation& sim, const Pos& pSource);
    bool PossibleRecruitTcyt(GrSimulation& sim, const Pos& pSource);
    bool PossibleRecruitTreg(GrSimulation& sim, const Pos& pSource);

};

inline bool RecruitmentProb::intCompareGTEQ(const double param1, const double param2)
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

inline bool RecruitmentProb::MacThresholdRecNew(const GrGrid &grid, const Pos &pSource, double &rThreshold)
{
    const double VmaxCCL2 = (1.0)/(3.0); // Temporary Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1
    const double VmaxCCL5 = (1.0)/(3.0); // Temporary Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1
    const double VmaxTNF = (1.0)/(3.0); // Temporary Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1


    const Scalar ratioCCL5toCCL2 = _PARAM(PARAM_MAC_SEC_RATE_CCL5) / _PARAM(PARAM_MAC_SEC_RATE_CCL2);

    rThreshold = (((VmaxCCL5 * (grid.CCL5(pSource) - (_PARAM(PARAM_TGAM_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCCL5toCCL2)))/((grid.CCL5(pSource) - (_PARAM(PARAM_TGAM_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCCL5toCCL2)) + (ratioCCL5toCCL2 * _PARAM(PARAM_TGAM_HALF_SAT_CHEMOKINE)))) +
                         ((VmaxTNF * (grid.TNF(pSource) - _PARAM(PARAM_TGAM_THRESHOLD_TNF_RECRUITMENT)))/((grid.TNF(pSource) - _PARAM(PARAM_TGAM_THRESHOLD_TNF_RECRUITMENT)) + _PARAM(PARAM_TGAM_HALF_SAT_TNF))) +
                                ((VmaxCCL2 * (grid.CCL2(pSource) - (_PARAM(PARAM_TGAM_THRESHOLD_CHEMOKINE_RECRUITMENT))))/((grid.CCL2(pSource) - (_PARAM(PARAM_TGAM_THRESHOLD_CHEMOKINE_RECRUITMENT))) + (_PARAM(PARAM_TGAM_HALF_SAT_CHEMOKINE)))));

//    rThreshold = (((VmaxCCL2 * grid.CCL2(pSource))/(grid.CCL2(pSource) + _PARAM(PARAM_MAC_HALF_SAT_CHEMOKINE))) +
//                         ((VmaxCCL5 * grid.CCL5(pSource))/(grid.CCL5(pSource) + (ratioCCL5toCCL2 * _PARAM(PARAM_MAC_HALF_SAT_CHEMOKINE)))) +
//                         ((VmaxTNF * grid.TNF(pSource))/(grid.TNF(pSource) + _PARAM(PARAM_MAC_HALF_SAT_TNF))));


    return (intCompareGTEQ(grid.CCL5(pSource), (_PARAM(PARAM_TCYT_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCCL5toCCL2)) && intCompareGTEQ(grid.TNF(pSource), _PARAM(PARAM_TCYT_THRESHOLD_TNF_RECRUITMENT)) && intCompareGTEQ(grid.CCL2(pSource), (_PARAM(PARAM_TCYT_THRESHOLD_CHEMOKINE_RECRUITMENT))));

//    return intCompareGTEQ(rThreshold, _PARAM(PARAM_MAC_THRESHOLD_RECRUITMENT));

}

inline bool RecruitmentProb::TgamThresholdRecNew(const GrGrid &grid, const Pos &pSource, double &rThreshold)
{
    const double VmaxCCL2 = (1.0)/(4.0); // Temporary Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1
    const double VmaxCCL5 = (1.0)/(4.0); // Temporary Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1
    const double VmaxTNF = (1.0)/(4.0); // Temporary Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1
    const double VmaxCXCL9 = (1.0)/(4.0); // Temporary Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1


    const Scalar ratioCCL5toCCL2 = _PARAM(PARAM_MAC_SEC_RATE_CCL5) / _PARAM(PARAM_MAC_SEC_RATE_CCL2);
    const Scalar ratioCXCL9toCCL2 = _PARAM(PARAM_MAC_SEC_RATE_CXCL9) / _PARAM(PARAM_MAC_SEC_RATE_CCL2);

//    rThreshold = (((VmaxCCL2 * grid.CCL2(pSource))/(grid.CCL2(pSource) + _PARAM(PARAM_TGAM_HALF_SAT_CHEMOKINE))) +
//                         ((VmaxCCL5 * grid.CCL5(pSource))/(grid.CCL5(pSource) + (ratioCCL5toCCL2 * _PARAM(PARAM_TGAM_HALF_SAT_CHEMOKINE)))) +
//                         ((VmaxTNF * grid.TNF(pSource))/(grid.TNF(pSource) + _PARAM(PARAM_TGAM_HALF_SAT_TNF))) +
//                          ((VmaxCXCL9 * grid.CXCL9(pSource))/(grid.CXCL9(pSource) + (ratioCXCL9toCCL2 * _PARAM(PARAM_TGAM_HALF_SAT_CHEMOKINE)))));


    rThreshold = (((VmaxCCL5 * (grid.CCL5(pSource) - (_PARAM(PARAM_TGAM_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCCL5toCCL2)))/((grid.CCL5(pSource) - (_PARAM(PARAM_TGAM_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCCL5toCCL2)) + (ratioCCL5toCCL2 * _PARAM(PARAM_TGAM_HALF_SAT_CHEMOKINE)))) +
                         ((VmaxTNF * (grid.TNF(pSource) - _PARAM(PARAM_TGAM_THRESHOLD_TNF_RECRUITMENT)))/((grid.TNF(pSource) - _PARAM(PARAM_TGAM_THRESHOLD_TNF_RECRUITMENT)) + _PARAM(PARAM_TGAM_HALF_SAT_TNF))) +
                            ((VmaxCXCL9 * (grid.CXCL9(pSource) - (_PARAM(PARAM_TGAM_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCXCL9toCCL2)))/((grid.CXCL9(pSource) - (_PARAM(PARAM_TGAM_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCXCL9toCCL2)) + (ratioCXCL9toCCL2 * _PARAM(PARAM_TGAM_HALF_SAT_CHEMOKINE)))) +
                                ((VmaxCCL2 * (grid.CCL2(pSource) - (_PARAM(PARAM_TGAM_THRESHOLD_CHEMOKINE_RECRUITMENT))))/((grid.CCL2(pSource) - (_PARAM(PARAM_TGAM_THRESHOLD_CHEMOKINE_RECRUITMENT))) + (_PARAM(PARAM_TGAM_HALF_SAT_CHEMOKINE)))));

//    std::cout << "Tgam Recruitment" << std::endl;
//    std::cout << pSource << "   " << rThreshold << std::endl;

    return (intCompareGTEQ(grid.CCL5(pSource), (_PARAM(PARAM_TCYT_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCCL5toCCL2)) && intCompareGTEQ(grid.TNF(pSource), _PARAM(PARAM_TCYT_THRESHOLD_TNF_RECRUITMENT)) && intCompareGTEQ(grid.CXCL9(pSource), (_PARAM(PARAM_TCYT_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCXCL9toCCL2)) && intCompareGTEQ(grid.CCL2(pSource), (_PARAM(PARAM_TCYT_THRESHOLD_CHEMOKINE_RECRUITMENT))));

//    return intCompareGTEQ(rThreshold, _PARAM(PARAM_TGAM_THRESHOLD_RECRUITMENT));
}

inline bool RecruitmentProb::TcytThresholdRecNew(const GrGrid &grid, const Pos &pSource, double &rThreshold)
{
    const double VmaxCCL5 = (1.0)/(3.0); // Temporary Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1
    const double VmaxTNF = (1.0)/(3.0); // Temporary Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1
    const double VmaxCXCL9 = (1.0)/(3.0); // Temporary Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1


    const Scalar ratioCCL5toCCL2 = _PARAM(PARAM_MAC_SEC_RATE_CCL5) / _PARAM(PARAM_MAC_SEC_RATE_CCL2);
    const Scalar ratioCXCL9toCCL2 = _PARAM(PARAM_MAC_SEC_RATE_CXCL9) / _PARAM(PARAM_MAC_SEC_RATE_CCL2);

//    rThreshold = (((VmaxCCL5 * grid.CCL5(pSource))/(grid.CCL5(pSource) + (ratioCCL5toCCL2 * _PARAM(PARAM_TCYT_HALF_SAT_CHEMOKINE)))) +
//                         ((VmaxTNF * grid.TNF(pSource))/(grid.TNF(pSource) + _PARAM(PARAM_TCYT_HALF_SAT_TNF))) +
//                          ((VmaxCXCL9 * grid.CXCL9(pSource))/(grid.CXCL9(pSource) + (ratioCXCL9toCCL2 * _PARAM(PARAM_TCYT_HALF_SAT_CHEMOKINE)))));


    rThreshold = (((VmaxCCL5 * (grid.CCL5(pSource) - (_PARAM(PARAM_TCYT_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCCL5toCCL2)))/((grid.CCL5(pSource) - (_PARAM(PARAM_TCYT_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCCL5toCCL2)) + (ratioCCL5toCCL2 * _PARAM(PARAM_TCYT_HALF_SAT_CHEMOKINE)))) +
                         ((VmaxTNF * (grid.TNF(pSource) - _PARAM(PARAM_TCYT_THRESHOLD_TNF_RECRUITMENT)))/((grid.TNF(pSource) - _PARAM(PARAM_TCYT_THRESHOLD_TNF_RECRUITMENT)) + _PARAM(PARAM_TCYT_HALF_SAT_TNF))) +
                            ((VmaxCXCL9 * (grid.CXCL9(pSource) - (_PARAM(PARAM_TCYT_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCXCL9toCCL2)))/((grid.CXCL9(pSource) - (_PARAM(PARAM_TCYT_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCXCL9toCCL2)) + (ratioCXCL9toCCL2 * _PARAM(PARAM_TCYT_HALF_SAT_CHEMOKINE)))));


//    std::cout << "Tcyt Recruitment" << std::endl;
//    std::cout << pSource << "   " << rThreshold << std::endl;

    return (intCompareGTEQ(grid.CCL5(pSource), (_PARAM(PARAM_TCYT_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCCL5toCCL2)) && intCompareGTEQ(grid.TNF(pSource), _PARAM(PARAM_TCYT_THRESHOLD_TNF_RECRUITMENT)) && intCompareGTEQ(grid.CXCL9(pSource), (_PARAM(PARAM_TCYT_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCXCL9toCCL2)));

//    return intCompareGTEQ(rThreshold, _PARAM(PARAM_TCYT_THRESHOLD_RECRUITMENT));
}

inline bool RecruitmentProb::TregThresholdRecNew(const GrGrid &grid, const Pos &pSource, double &rThreshold)
{
    const double VmaxCCL5 = (1.0)/(4.0); // Temporary Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1
    const double VmaxTNF = (3.0)/(4.0); // Temporary Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1

    const Scalar ratioCCL5toCCL2 = _PARAM(PARAM_MAC_SEC_RATE_CCL5) / _PARAM(PARAM_MAC_SEC_RATE_CCL2);

    rThreshold = (((VmaxCCL5 * (grid.CCL5(pSource) - (_PARAM(PARAM_TREG_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCCL5toCCL2)))/((grid.CCL5(pSource) - (_PARAM(PARAM_TREG_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCCL5toCCL2)) + (ratioCCL5toCCL2 * _PARAM(PARAM_TREG_HALF_SAT_CHEMOKINE)))) +
                         ((VmaxTNF * (grid.TNF(pSource) - _PARAM(PARAM_TREG_THRESHOLD_TNF_RECRUITMENT)))/((grid.TNF(pSource) - _PARAM(PARAM_TREG_THRESHOLD_TNF_RECRUITMENT)) + _PARAM(PARAM_TREG_HALF_SAT_TNF))));

//    std::cout << "Treg Recruitment" << std::endl;
//    std::cout << pSource << "   " << rThreshold << std::endl;

    return (intCompareGTEQ(grid.CCL5(pSource), (_PARAM(PARAM_TREG_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCCL5toCCL2)) && intCompareGTEQ(grid.TNF(pSource), _PARAM(PARAM_TREG_THRESHOLD_TNF_RECRUITMENT)));

//    return intCompareGTEQ(rThreshold, _PARAM(PARAM_TREG_THRESHOLD_RECRUITMENT));
}

inline bool RecruitmentProb::MacThresholdRecNew(const GrGrid &grid, const Pos &pSource)
{
    double threshold;
    return MacThresholdRecNew(grid, pSource, threshold);
}

inline bool RecruitmentProb::TgamThresholdRecNew(const GrGrid &grid, const Pos &pSource)
{
    double threshold;
    return TgamThresholdRecNew(grid, pSource, threshold);
}

inline bool RecruitmentProb::TcytThresholdRecNew(const GrGrid &grid, const Pos &pSource)
{
    double threshold;
    return TcytThresholdRecNew(grid, pSource, threshold);
}

inline bool RecruitmentProb::TregThresholdRecNew(const GrGrid &grid, const Pos &pSource)
{
    double threshold;
    return TregThresholdRecNew(grid, pSource, threshold);
}

inline RecruitmentMethod RecruitmentProb::getMethod() const
{
	return RECR_PROB;
}

inline void RecruitmentProb::serialize(std::ostream&) const
{
}

inline void RecruitmentProb::deserialize(std::istream&)
{
}
#endif /* RECRUITMENTPROB_H_ */
