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
    void recruit(GrSimulation& sim, int time);
    void recruit(GrSimulation& sim, const Pos& pSource, size_t cellRecNum);

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

    bool PossibleRecruitMac(GrSimulation& sim, const Pos& pSource, size_t cellRecNum);
    bool PossibleRecruitTgam(GrSimulation& sim, const Pos& pSource, size_t cellRecNum);
    bool PossibleRecruitTcyt(GrSimulation& sim, const Pos& pSource, size_t cellRecNum);
    bool PossibleRecruitTreg(GrSimulation& sim, const Pos& pSource, size_t cellRecNum);

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

inline bool RecruitmentProb::MacThresholdRecNew(const GrGrid &grid, const Pos &pSource, double &rThreshold)
{
    const double VmaxCCL2 = (1.0)/(3.0); // Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1
    const double VmaxCCL5 = (1.0)/(3.0); // Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1
    const double VmaxTNF = (1.0)/(3.0); // Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1

    const Scalar ratioCCL5toCCL2 = _PARAM(PARAM_MAC_SEC_RATE_CCL5) / _PARAM(PARAM_MAC_SEC_RATE_CCL2);

    bool thresholdCCL2 = intCompareGTEQ(grid.CCL2(pSource), (_PARAM(PARAM_MAC_THRESHOLD_CHEMOKINE_RECRUITMENT)));
    bool thresholdCCL5 = intCompareGTEQ(grid.CCL5(pSource), (_PARAM(PARAM_MAC_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCCL5toCCL2));
    bool thresholdTNF = intCompareGTEQ(grid.TNF(pSource), _PARAM(PARAM_MAC_THRESHOLD_TNF_RECRUITMENT));

    double rThresholdCCL5 = 0.0, rThresholdTNF = 0.0, rThresholdCCL2 = 0.0;

    if (thresholdCCL2)
        rThresholdCCL2 = ((VmaxCCL2 * (grid.CCL2(pSource) - (_PARAM(PARAM_MAC_THRESHOLD_CHEMOKINE_RECRUITMENT))))/((grid.CCL2(pSource) - (_PARAM(PARAM_MAC_THRESHOLD_CHEMOKINE_RECRUITMENT))) + (_PARAM(PARAM_MAC_HALF_SAT_CHEMOKINE))));

    if (thresholdCCL5)
        rThresholdCCL5 = ((VmaxCCL5 * (grid.CCL5(pSource) - (_PARAM(PARAM_MAC_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCCL5toCCL2)))/((grid.CCL5(pSource) - (_PARAM(PARAM_MAC_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCCL5toCCL2)) + (ratioCCL5toCCL2 * _PARAM(PARAM_MAC_HALF_SAT_CHEMOKINE))));

    if (thresholdTNF)
        rThresholdTNF = ((VmaxTNF * (grid.TNF(pSource) - _PARAM(PARAM_MAC_THRESHOLD_TNF_RECRUITMENT)))/((grid.TNF(pSource) - _PARAM(PARAM_MAC_THRESHOLD_TNF_RECRUITMENT)) + _PARAM(PARAM_MAC_HALF_SAT_TNF)));

    rThreshold = rThresholdCCL2 + rThresholdCCL5 + rThresholdTNF;

//    rThreshold = (thresholdCCL5 * (((VmaxCCL5 * (grid.CCL5(pSource) - (_PARAM(PARAM_MAC_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCCL5toCCL2)))/((grid.CCL5(pSource) - (_PARAM(PARAM_MAC_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCCL5toCCL2)) + (ratioCCL5toCCL2 * _PARAM(PARAM_MAC_HALF_SAT_CHEMOKINE))))) +
//                         (thresholdTNF * ((VmaxTNF * (grid.TNF(pSource) - _PARAM(PARAM_MAC_THRESHOLD_TNF_RECRUITMENT)))/((grid.TNF(pSource) - _PARAM(PARAM_MAC_THRESHOLD_TNF_RECRUITMENT)) + _PARAM(PARAM_MAC_HALF_SAT_TNF)))) +
//                                (thresholdCCL2 * ((VmaxCCL2 * (grid.CCL2(pSource) - (_PARAM(PARAM_MAC_THRESHOLD_CHEMOKINE_RECRUITMENT))))/((grid.CCL2(pSource) - (_PARAM(PARAM_MAC_THRESHOLD_CHEMOKINE_RECRUITMENT))) + (_PARAM(PARAM_MAC_HALF_SAT_CHEMOKINE))))));

    return (thresholdCCL5 || thresholdTNF || thresholdCCL2);
}

inline bool RecruitmentProb::TgamThresholdRecNew(const GrGrid &grid, const Pos &pSource, double &rThreshold)
{
    const double VmaxCCL2 = (1.0)/(4.0); // Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1
    const double VmaxCCL5 = (1.0)/(4.0); // Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1
    const double VmaxTNF = (1.0)/(4.0); // Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1
    const double VmaxCXCL9 = (1.0)/(4.0); // Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1

    const Scalar ratioCCL5toCCL2 = _PARAM(PARAM_MAC_SEC_RATE_CCL5) / _PARAM(PARAM_MAC_SEC_RATE_CCL2);
    const Scalar ratioCXCL9toCCL2 = _PARAM(PARAM_MAC_SEC_RATE_CXCL9) / _PARAM(PARAM_MAC_SEC_RATE_CCL2);

    bool thresholdCCL5 = intCompareGTEQ(grid.CCL5(pSource), (_PARAM(PARAM_TGAM_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCCL5toCCL2));
    bool thresholdTNF = intCompareGTEQ(grid.TNF(pSource), _PARAM(PARAM_TGAM_THRESHOLD_TNF_RECRUITMENT));
    bool thresholdCXCL9 = intCompareGTEQ(grid.CXCL9(pSource), (_PARAM(PARAM_TGAM_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCXCL9toCCL2));
    bool thresholdCCL2 = intCompareGTEQ(grid.CCL2(pSource), (_PARAM(PARAM_TGAM_THRESHOLD_CHEMOKINE_RECRUITMENT)));

    double rThresholdCCL5 = 0.0, rThresholdTNF = 0.0, rThresholdCXCL9 = 0.0, rThresholdCCL2 = 0.0;

    if (thresholdCCL5)
        rThresholdCCL5 = ((VmaxCCL5 * (grid.CCL5(pSource) - (_PARAM(PARAM_TGAM_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCCL5toCCL2)))/((grid.CCL5(pSource) - (_PARAM(PARAM_TGAM_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCCL5toCCL2)) + (ratioCCL5toCCL2 * _PARAM(PARAM_TGAM_HALF_SAT_CHEMOKINE))));

    if (thresholdTNF)
        rThresholdTNF = ((VmaxTNF * (grid.TNF(pSource) - _PARAM(PARAM_TGAM_THRESHOLD_TNF_RECRUITMENT)))/((grid.TNF(pSource) - _PARAM(PARAM_TGAM_THRESHOLD_TNF_RECRUITMENT)) + _PARAM(PARAM_TGAM_HALF_SAT_TNF)));

    if (thresholdCXCL9)
        rThresholdCXCL9 = ((VmaxCXCL9 * (grid.CXCL9(pSource) - (_PARAM(PARAM_TGAM_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCXCL9toCCL2)))/((grid.CXCL9(pSource) - (_PARAM(PARAM_TGAM_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCXCL9toCCL2)) + (ratioCXCL9toCCL2 * _PARAM(PARAM_TGAM_HALF_SAT_CHEMOKINE))));

    if (thresholdCCL2)
        rThresholdCCL2 = ((VmaxCCL2 * (grid.CCL2(pSource) - (_PARAM(PARAM_TGAM_THRESHOLD_CHEMOKINE_RECRUITMENT))))/((grid.CCL2(pSource) - (_PARAM(PARAM_TGAM_THRESHOLD_CHEMOKINE_RECRUITMENT))) + (_PARAM(PARAM_TGAM_HALF_SAT_CHEMOKINE))));

    rThreshold = rThresholdCCL5 + rThresholdTNF + rThresholdCXCL9 + rThresholdCCL2;

//    rThreshold = (thresholdCCL5 * ((VmaxCCL5 * (grid.CCL5(pSource) - (_PARAM(PARAM_TGAM_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCCL5toCCL2)))/((grid.CCL5(pSource) - (_PARAM(PARAM_TGAM_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCCL5toCCL2)) + (ratioCCL5toCCL2 * _PARAM(PARAM_TGAM_HALF_SAT_CHEMOKINE))))) +
//                         (thresholdTNF * ((VmaxTNF * (grid.TNF(pSource) - _PARAM(PARAM_TGAM_THRESHOLD_TNF_RECRUITMENT)))/((grid.TNF(pSource) - _PARAM(PARAM_TGAM_THRESHOLD_TNF_RECRUITMENT)) + _PARAM(PARAM_TGAM_HALF_SAT_TNF)))) +
//                            (thresholdCXCL9 * ((VmaxCXCL9 * (grid.CXCL9(pSource) - (_PARAM(PARAM_TGAM_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCXCL9toCCL2)))/((grid.CXCL9(pSource) - (_PARAM(PARAM_TGAM_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCXCL9toCCL2)) + (ratioCXCL9toCCL2 * _PARAM(PARAM_TGAM_HALF_SAT_CHEMOKINE))))) +
//                                (thresholdCCL2 * ((VmaxCCL2 * (grid.CCL2(pSource) - (_PARAM(PARAM_TGAM_THRESHOLD_CHEMOKINE_RECRUITMENT))))/((grid.CCL2(pSource) - (_PARAM(PARAM_TGAM_THRESHOLD_CHEMOKINE_RECRUITMENT))) + (_PARAM(PARAM_TGAM_HALF_SAT_CHEMOKINE)))));

    return (thresholdCCL5 || thresholdTNF || thresholdCXCL9 || thresholdCCL2);
}

inline bool RecruitmentProb::TcytThresholdRecNew(const GrGrid &grid, const Pos &pSource, double &rThreshold)
{
    const double VmaxCCL5 = (1.0)/(3.0); // Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1
    const double VmaxTNF = (1.0)/(3.0); // Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1
    const double VmaxCXCL9 = (1.0)/(3.0); // Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1

    const Scalar ratioCCL5toCCL2 = _PARAM(PARAM_MAC_SEC_RATE_CCL5) / _PARAM(PARAM_MAC_SEC_RATE_CCL2);
    const Scalar ratioCXCL9toCCL2 = _PARAM(PARAM_MAC_SEC_RATE_CXCL9) / _PARAM(PARAM_MAC_SEC_RATE_CCL2);

    bool thresholdCCL5 = intCompareGTEQ(grid.CCL5(pSource), (_PARAM(PARAM_TCYT_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCCL5toCCL2));
    bool thresholdTNF = intCompareGTEQ(grid.TNF(pSource), _PARAM(PARAM_TCYT_THRESHOLD_TNF_RECRUITMENT));
    bool thresholdCXCL9 = intCompareGTEQ(grid.CXCL9(pSource), (_PARAM(PARAM_TCYT_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCXCL9toCCL2));

    double rThresholdCCL5 = 0.0, rThresholdTNF = 0.0, rThresholdCXCL9 = 0.0;

    if (thresholdCCL5)
        rThresholdCCL5 = ((VmaxCCL5 * (grid.CCL5(pSource) - (_PARAM(PARAM_TCYT_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCCL5toCCL2)))/((grid.CCL5(pSource) - (_PARAM(PARAM_TCYT_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCCL5toCCL2)) + (ratioCCL5toCCL2 * _PARAM(PARAM_TCYT_HALF_SAT_CHEMOKINE))));

    if (thresholdTNF)
        rThresholdTNF = ((VmaxTNF * (grid.TNF(pSource) - _PARAM(PARAM_TCYT_THRESHOLD_TNF_RECRUITMENT)))/((grid.TNF(pSource) - _PARAM(PARAM_TCYT_THRESHOLD_TNF_RECRUITMENT)) + _PARAM(PARAM_TCYT_HALF_SAT_TNF)));

    if (thresholdCXCL9)
        rThresholdCXCL9 = ((VmaxCXCL9 * (grid.CXCL9(pSource) - (_PARAM(PARAM_TCYT_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCXCL9toCCL2)))/((grid.CXCL9(pSource) - (_PARAM(PARAM_TCYT_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCXCL9toCCL2)) + (ratioCXCL9toCCL2 * _PARAM(PARAM_TCYT_HALF_SAT_CHEMOKINE))));

    rThreshold = rThresholdCCL5 + rThresholdTNF + rThresholdCXCL9;

//    rThreshold = (thresholdCCL5 * ((VmaxCCL5 * (grid.CCL5(pSource) - (_PARAM(PARAM_TCYT_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCCL5toCCL2)))/((grid.CCL5(pSource) - (_PARAM(PARAM_TCYT_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCCL5toCCL2)) + (ratioCCL5toCCL2 * _PARAM(PARAM_TCYT_HALF_SAT_CHEMOKINE))))) +
//                         (thresholdTNF * ((VmaxTNF * (grid.TNF(pSource) - _PARAM(PARAM_TCYT_THRESHOLD_TNF_RECRUITMENT)))/((grid.TNF(pSource) - _PARAM(PARAM_TCYT_THRESHOLD_TNF_RECRUITMENT)) + _PARAM(PARAM_TCYT_HALF_SAT_TNF)))) +
//                            (thresholdCXCL9 * ((VmaxCXCL9 * (grid.CXCL9(pSource) - (_PARAM(PARAM_TCYT_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCXCL9toCCL2)))/((grid.CXCL9(pSource) - (_PARAM(PARAM_TCYT_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCXCL9toCCL2)) + (ratioCXCL9toCCL2 * _PARAM(PARAM_TCYT_HALF_SAT_CHEMOKINE)))));

    return (thresholdCCL5 || thresholdTNF || thresholdCXCL9);
}

inline bool RecruitmentProb::TregThresholdRecNew(const GrGrid &grid, const Pos &pSource, double &rThreshold)
{
    const double VmaxCCL5 = (1.0)/(4.0); // Temporary Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1
    const double VmaxTNF = (3.0)/(4.0); // Temporary Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1

    const Scalar ratioCCL5toCCL2 = _PARAM(PARAM_MAC_SEC_RATE_CCL5) / _PARAM(PARAM_MAC_SEC_RATE_CCL2);

    bool thresholdCCL5 = intCompareGTEQ(grid.CCL5(pSource), (_PARAM(PARAM_TREG_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCCL5toCCL2));
    bool thresholdTNF = intCompareGTEQ(grid.TNF(pSource), _PARAM(PARAM_TREG_THRESHOLD_TNF_RECRUITMENT));

    double rThresholdCCL5 = 0.0, rThresholdTNF = 0.0;

    if (thresholdCCL5)
        rThresholdCCL5 = ((VmaxCCL5 * (grid.CCL5(pSource) - (_PARAM(PARAM_TREG_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCCL5toCCL2)))/((grid.CCL5(pSource) - (_PARAM(PARAM_TREG_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCCL5toCCL2)) + (ratioCCL5toCCL2 * _PARAM(PARAM_TREG_HALF_SAT_CHEMOKINE))));

    if (thresholdTNF)
        rThresholdTNF = ((VmaxTNF * (grid.TNF(pSource) - _PARAM(PARAM_TREG_THRESHOLD_TNF_RECRUITMENT)))/((grid.TNF(pSource) - _PARAM(PARAM_TREG_THRESHOLD_TNF_RECRUITMENT)) + _PARAM(PARAM_TREG_HALF_SAT_TNF)));

    rThreshold = rThresholdCCL5 + rThresholdTNF;

//    rThreshold = (thresholdCCL5 * ((VmaxCCL5 * (grid.CCL5(pSource) - (_PARAM(PARAM_TREG_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCCL5toCCL2)))/((grid.CCL5(pSource) - (_PARAM(PARAM_TREG_THRESHOLD_CHEMOKINE_RECRUITMENT) * ratioCCL5toCCL2)) + (ratioCCL5toCCL2 * _PARAM(PARAM_TREG_HALF_SAT_CHEMOKINE))))) +
//                         (thresholdTNF * ((VmaxTNF * (grid.TNF(pSource) - _PARAM(PARAM_TREG_THRESHOLD_TNF_RECRUITMENT)))/((grid.TNF(pSource) - _PARAM(PARAM_TREG_THRESHOLD_TNF_RECRUITMENT)) + _PARAM(PARAM_TREG_HALF_SAT_TNF))));

    return (thresholdCCL5 || thresholdTNF);
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
