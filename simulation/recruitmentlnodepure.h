/*
 * RecruitmentLnODEPure.h
 *
 *  Created on: Jul 4, 2010
 *      Author: mohammed
 */

#ifndef RECRUITMENTLNODEPURE_H_
#define RECRUITMENTLNODEPURE_H_

#include "recruitmentlnode.h"

class RecruitmentLnODEPure : public RecruitmentLnODE
{
private:
  static const int _idxTH0 = 0;
  static const int _idxTH1 = 1;
  static const int _idxT80 = 2;
  static const int _idxT8 = 3;
  static const int _idxPrecursorCTL = 4;
  static const int _idxEffectorTH0 = 5;
  static const int _idxEffectorT80 = 6;

  virtual void solveODE(const int time, const Stats& statsPrevious, Stats& stats);

public:
  RecruitmentLnODEPure();
  virtual ~RecruitmentLnODEPure();

  RecruitmentMethod getMethod() const;
  /**
  * @copydoc GrSimulation::serialize
  */
  template<class Archive>
  void serialize(Archive& ar, const unsigned int version);

  /*virtual*/
  RecruitmentBase* clone() const
  {
    return new RecruitmentLnODEPure(*this);
  }

};

inline RecruitmentMethod RecruitmentLnODEPure::getMethod() const
{
  return RECR_LN_ODE_PURE;
}

template<class Archive>
inline void RecruitmentLnODEPure::serialize(Archive& ar, const unsigned int /*version*/)
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(RecruitmentLnODE);
}

BOOST_CLASS_EXPORT_KEY(RecruitmentLnODEPure)
#endif /* RECRUITMENTLNODEPURE_H_ */
