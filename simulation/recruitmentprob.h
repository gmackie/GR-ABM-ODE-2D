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
  /**
  * @copydoc GrSimulation::serialize
  */
  template<class Archive>
  void serialize(Archive& ar, const unsigned int version);


  static bool MacThresholdRecNew(const GrGrid& grid, const Pos& pSource);
  static bool TgamThresholdRecNew(const GrGrid& grid, const Pos& pSource);
  static bool TcytThresholdRecNew(const GrGrid& grid, const Pos& pSource);
  static bool TregThresholdRecNew(const GrGrid& grid, const Pos& pSource);

  /*virtual*/
  RecruitmentBase* clone() const
  {
    return new RecruitmentProb(*this);
  }
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

inline bool RecruitmentProb::MacThresholdRecNew(const GrGrid &grid, const Pos &pSource, double &rThreshold)
{
  const double VmaxCCL2 = (3.0)/(8.0); // Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1
  const double VmaxCCL5 = (3.0)/(8.0); // Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1
  const double VmaxTNF = (1.0)/(4.0); // Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1

  const Scalar ratioCCL5toCCL2 = _PARAM(Mac_dCCL5) / _PARAM(Mac_dCCL2);

  bool thresholdCCL2 = intCompareGTEQ(grid.CCL2(pSource), (_PARAM(Mac_thresholdRecChemokine)));
  bool thresholdCCL5 = intCompareGTEQ(grid.CCL5(pSource), (_PARAM(Mac_thresholdRecChemokine) * ratioCCL5toCCL2));
  bool thresholdTNF = intCompareGTEQ(grid.TNF(pSource), _PARAM(Mac_thresholdRecTNF));

  double rThresholdCCL5 = 0.0, rThresholdTNF = 0.0, rThresholdCCL2 = 0.0;

  if (thresholdCCL2)
    rThresholdCCL2 = ((VmaxCCL2 * (grid.CCL2(pSource) - (_PARAM(Mac_thresholdRecChemokine))))/((grid.CCL2(pSource) - (_PARAM(Mac_thresholdRecChemokine))) + (_PARAM(Mac_recruitmentHalfSatChemokine))));

  if (thresholdCCL5)
    rThresholdCCL5 = ((VmaxCCL5 * (grid.CCL5(pSource) - (_PARAM(Mac_thresholdRecChemokine) * ratioCCL5toCCL2)))/((grid.CCL5(pSource) - (_PARAM(Mac_thresholdRecChemokine) * ratioCCL5toCCL2)) + (ratioCCL5toCCL2 * _PARAM(Mac_recruitmentHalfSatChemokine))));

  if (thresholdTNF)
    rThresholdTNF = ((VmaxTNF * (grid.TNF(pSource) - _PARAM(Mac_thresholdRecTNF)))/((grid.TNF(pSource) - _PARAM(Mac_thresholdRecTNF)) + _PARAM(Mac_recruitmentHalfSatTNF)));

  rThreshold = rThresholdCCL2 + rThresholdCCL5 + rThresholdTNF;

//    rThreshold = (thresholdCCL5 * (((VmaxCCL5 * (grid.CCL5(pSource) - (_PARAM(Mac_thresholdRecChemokine) * ratioCCL5toCCL2)))/((grid.CCL5(pSource) - (_PARAM(Mac_thresholdRecChemokine) * ratioCCL5toCCL2)) + (ratioCCL5toCCL2 * _PARAM(Mac_recruitmentHalfSatChemokine))))) +
//                         (thresholdTNF * ((VmaxTNF * (grid.TNF(pSource) - _PARAM(Mac_thresholdRecTNF)))/((grid.TNF(pSource) - _PARAM(Mac_thresholdRecTNF)) + _PARAM(Mac_recruitmentHalfSatTNF)))) +
//                                (thresholdCCL2 * ((VmaxCCL2 * (grid.CCL2(pSource) - (_PARAM(Mac_thresholdRecChemokine))))/((grid.CCL2(pSource) - (_PARAM(Mac_thresholdRecChemokine))) + (_PARAM(Mac_recruitmentHalfSatChemokine))))));

  return (thresholdCCL5 || thresholdTNF || thresholdCCL2);
}

inline bool RecruitmentProb::TgamThresholdRecNew(const GrGrid &grid, const Pos &pSource, double &rThreshold)
{
  const double VmaxCCL2 = (2.0)/(9.0); // Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1
  const double VmaxCCL5 = (2.0)/(9.0); // Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1
  const double VmaxTNF = (1.0)/(3.0); // Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1
  const double VmaxCXCL9 = (2.0)/(9.0); // Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1

  const Scalar ratioCCL5toCCL2 = _PARAM(Mac_dCCL5) / _PARAM(Mac_dCCL2);
  const Scalar ratioCXCL9toCCL2 = _PARAM(Mac_dCXCL9) / _PARAM(Mac_dCCL2);

  bool thresholdCCL5 = intCompareGTEQ(grid.CCL5(pSource), (_PARAM(Tcell_Tgam_thresholdRecChemokine) * ratioCCL5toCCL2));
  bool thresholdTNF = intCompareGTEQ(grid.TNF(pSource), _PARAM(Tcell_Tgam_thresholdRecTNF));
  bool thresholdCXCL9 = intCompareGTEQ(grid.CXCL9(pSource), (_PARAM(Tcell_Tgam_thresholdRecChemokine) * ratioCXCL9toCCL2));
  bool thresholdCCL2 = intCompareGTEQ(grid.CCL2(pSource), (_PARAM(Tcell_Tgam_thresholdRecChemokine)));

  double rThresholdCCL5 = 0.0, rThresholdTNF = 0.0, rThresholdCXCL9 = 0.0, rThresholdCCL2 = 0.0;

  if (thresholdCCL5)
    rThresholdCCL5 = ((VmaxCCL5 * (grid.CCL5(pSource) - (_PARAM(Tcell_Tgam_thresholdRecChemokine) * ratioCCL5toCCL2)))/((grid.CCL5(pSource) - (_PARAM(Tcell_Tgam_thresholdRecChemokine) * ratioCCL5toCCL2)) + (ratioCCL5toCCL2 * _PARAM(Tcell_Tgam_recruitmentHalfSatChemokine))));

  if (thresholdTNF)
    rThresholdTNF = ((VmaxTNF * (grid.TNF(pSource) - _PARAM(Tcell_Tgam_thresholdRecTNF)))/((grid.TNF(pSource) - _PARAM(Tcell_Tgam_thresholdRecTNF)) + _PARAM(Tcell_Tgam_recruitmentHalfSatTNF)));

  if (thresholdCXCL9)
    rThresholdCXCL9 = ((VmaxCXCL9 * (grid.CXCL9(pSource) - (_PARAM(Tcell_Tgam_thresholdRecChemokine) * ratioCXCL9toCCL2)))/((grid.CXCL9(pSource) - (_PARAM(Tcell_Tgam_thresholdRecChemokine) * ratioCXCL9toCCL2)) + (ratioCXCL9toCCL2 * _PARAM(Tcell_Tgam_recruitmentHalfSatChemokine))));

  if (thresholdCCL2)
    rThresholdCCL2 = ((VmaxCCL2 * (grid.CCL2(pSource) - (_PARAM(Tcell_Tgam_thresholdRecChemokine))))/((grid.CCL2(pSource) - (_PARAM(Tcell_Tgam_thresholdRecChemokine))) + (_PARAM(Tcell_Tgam_recruitmentHalfSatChemokine))));

  rThreshold = rThresholdCCL5 + rThresholdTNF + rThresholdCXCL9 + rThresholdCCL2;

//    rThreshold = (thresholdCCL5 * ((VmaxCCL5 * (grid.CCL5(pSource) - (_PARAM(Tcell_Tgam_thresholdRecChemokine) * ratioCCL5toCCL2)))/((grid.CCL5(pSource) - (_PARAM(Tcell_Tgam_thresholdRecChemokine) * ratioCCL5toCCL2)) + (ratioCCL5toCCL2 * _PARAM(Tcell_Tgam_recruitmentHalfSatChemokine))))) +
//                         (thresholdTNF * ((VmaxTNF * (grid.TNF(pSource) - _PARAM(Tcell_Tgam_thresholdRecTNF)))/((grid.TNF(pSource) - _PARAM(Tcell_Tgam_thresholdRecTNF)) + _PARAM(Tcell_Tgam_recruitmentHalfSatTNF)))) +
//                            (thresholdCXCL9 * ((VmaxCXCL9 * (grid.CXCL9(pSource) - (_PARAM(Tcell_Tgam_thresholdRecChemokine) * ratioCXCL9toCCL2)))/((grid.CXCL9(pSource) - (_PARAM(Tcell_Tgam_thresholdRecChemokine) * ratioCXCL9toCCL2)) + (ratioCXCL9toCCL2 * _PARAM(Tcell_Tgam_recruitmentHalfSatChemokine))))) +
//                                (thresholdCCL2 * ((VmaxCCL2 * (grid.CCL2(pSource) - (_PARAM(Tcell_Tgam_thresholdRecChemokine))))/((grid.CCL2(pSource) - (_PARAM(Tcell_Tgam_thresholdRecChemokine))) + (_PARAM(Tcell_Tgam_recruitmentHalfSatChemokine)))));

  return (thresholdCCL5 || thresholdTNF || thresholdCXCL9 || thresholdCCL2);
}

inline bool RecruitmentProb::TcytThresholdRecNew(const GrGrid &grid, const Pos &pSource, double &rThreshold)
{
  const double VmaxCCL5 = (1.0)/(3.0); // Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1
  const double VmaxTNF = (1.0)/(3.0); // Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1
  const double VmaxCXCL9 = (1.0)/(3.0); // Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1

  const Scalar ratioCCL5toCCL2 = _PARAM(Mac_dCCL5) / _PARAM(Mac_dCCL2);
  const Scalar ratioCXCL9toCCL2 = _PARAM(Mac_dCXCL9) / _PARAM(Mac_dCCL2);

  bool thresholdCCL5 = intCompareGTEQ(grid.CCL5(pSource), (_PARAM(Tcell_Tcyt_thresholdRecChemokine) * ratioCCL5toCCL2));
  bool thresholdTNF = intCompareGTEQ(grid.TNF(pSource), _PARAM(Tcell_Tcyt_thresholdRecTNF));
  bool thresholdCXCL9 = intCompareGTEQ(grid.CXCL9(pSource), (_PARAM(Tcell_Tcyt_thresholdRecChemokine) * ratioCXCL9toCCL2));

  double rThresholdCCL5 = 0.0, rThresholdTNF = 0.0, rThresholdCXCL9 = 0.0;

  if (thresholdCCL5)
    rThresholdCCL5 = ((VmaxCCL5 * (grid.CCL5(pSource) - (_PARAM(Tcell_Tcyt_thresholdRecChemokine) * ratioCCL5toCCL2)))/((grid.CCL5(pSource) - (_PARAM(Tcell_Tcyt_thresholdRecChemokine) * ratioCCL5toCCL2)) + (ratioCCL5toCCL2 * _PARAM(Tcell_Tcyt_recruitmentHalfSatChemokine))));

  if (thresholdTNF)
    rThresholdTNF = ((VmaxTNF * (grid.TNF(pSource) - _PARAM(Tcell_Tcyt_thresholdRecTNF)))/((grid.TNF(pSource) - _PARAM(Tcell_Tcyt_thresholdRecTNF)) + _PARAM(Tcell_Tcyt_recruitmentHalfSatTNF)));

  if (thresholdCXCL9)
    rThresholdCXCL9 = ((VmaxCXCL9 * (grid.CXCL9(pSource) - (_PARAM(Tcell_Tcyt_thresholdRecChemokine) * ratioCXCL9toCCL2)))/((grid.CXCL9(pSource) - (_PARAM(Tcell_Tcyt_thresholdRecChemokine) * ratioCXCL9toCCL2)) + (ratioCXCL9toCCL2 * _PARAM(Tcell_Tcyt_recruitmentHalfSatChemokine))));

  rThreshold = rThresholdCCL5 + rThresholdTNF + rThresholdCXCL9;

//    rThreshold = (thresholdCCL5 * ((VmaxCCL5 * (grid.CCL5(pSource) - (_PARAM(Tcell_Tcyt_thresholdRecChemokine) * ratioCCL5toCCL2)))/((grid.CCL5(pSource) - (_PARAM(Tcell_Tcyt_thresholdRecChemokine) * ratioCCL5toCCL2)) + (ratioCCL5toCCL2 * _PARAM(Tcell_Tcyt_recruitmentHalfSatChemokine))))) +
//                         (thresholdTNF * ((VmaxTNF * (grid.TNF(pSource) - _PARAM(Tcell_Tcyt_thresholdRecTNF)))/((grid.TNF(pSource) - _PARAM(Tcell_Tcyt_thresholdRecTNF)) + _PARAM(Tcell_Tcyt_recruitmentHalfSatTNF)))) +
//                            (thresholdCXCL9 * ((VmaxCXCL9 * (grid.CXCL9(pSource) - (_PARAM(Tcell_Tcyt_thresholdRecChemokine) * ratioCXCL9toCCL2)))/((grid.CXCL9(pSource) - (_PARAM(Tcell_Tcyt_thresholdRecChemokine) * ratioCXCL9toCCL2)) + (ratioCXCL9toCCL2 * _PARAM(Tcell_Tcyt_recruitmentHalfSatChemokine)))));

  return (thresholdCCL5 || thresholdTNF || thresholdCXCL9);
}

inline bool RecruitmentProb::TregThresholdRecNew(const GrGrid &grid, const Pos &pSource, double &rThreshold)
{
  const double VmaxCCL5 = (1.0)/(4.0); // Temporary Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1
  const double VmaxTNF = (3.0)/(4.0); // Temporary Vmax is set by the number of species recruitment is based off so all values scale between 0 and 1

  const Scalar ratioCCL5toCCL2 = _PARAM(Mac_dCCL5) / _PARAM(Mac_dCCL2);

  bool thresholdCCL5 = intCompareGTEQ(grid.CCL5(pSource), (_PARAM(Tcell_Treg_thresholdRecChemokine) * ratioCCL5toCCL2));
  bool thresholdTNF = intCompareGTEQ(grid.TNF(pSource), _PARAM(Tcell_Treg_thresholdRecTNF));

  double rThresholdCCL5 = 0.0, rThresholdTNF = 0.0;

  if (thresholdCCL5)
    rThresholdCCL5 = ((VmaxCCL5 * (grid.CCL5(pSource) - (_PARAM(Tcell_Treg_thresholdRecChemokine) * ratioCCL5toCCL2)))/((grid.CCL5(pSource) - (_PARAM(Tcell_Treg_thresholdRecChemokine) * ratioCCL5toCCL2)) + (ratioCCL5toCCL2 * _PARAM(Tcell_Treg_recruitmentHalfSatChemokine))));

  if (thresholdTNF)
    rThresholdTNF = ((VmaxTNF * (grid.TNF(pSource) - _PARAM(Tcell_Treg_thresholdRecTNF)))/((grid.TNF(pSource) - _PARAM(Tcell_Treg_thresholdRecTNF)) + _PARAM(Tcell_Treg_recruitmentHalfSatTNF)));

  rThreshold = rThresholdCCL5 + rThresholdTNF;

//    rThreshold = (thresholdCCL5 * ((VmaxCCL5 * (grid.CCL5(pSource) - (_PARAM(Tcell_Treg_thresholdRecChemokine) * ratioCCL5toCCL2)))/((grid.CCL5(pSource) - (_PARAM(Tcell_Treg_thresholdRecChemokine) * ratioCCL5toCCL2)) + (ratioCCL5toCCL2 * _PARAM(Tcell_Treg_recruitmentHalfSatChemokine))))) +
//                         (thresholdTNF * ((VmaxTNF * (grid.TNF(pSource) - _PARAM(Tcell_Treg_thresholdRecTNF)))/((grid.TNF(pSource) - _PARAM(Tcell_Treg_thresholdRecTNF)) + _PARAM(Tcell_Treg_recruitmentHalfSatTNF))));

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

template<class Archive>
inline void RecruitmentProb::serialize(Archive& ar, const unsigned int /*version*/)
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(RecruitmentBase);
}

BOOST_CLASS_EXPORT_KEY(RecruitmentProb)
#endif /* RECRUITMENTPROB_H_ */
