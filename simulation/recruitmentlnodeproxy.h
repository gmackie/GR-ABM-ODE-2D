/*
 * recruitmentlnodeproxy.h
 *
 *  Created on: Mar 2, 2012
 *      Author: pwolberg
 */

#include "recruitmentlnode.h"

#ifndef RECRUITMENTLNODEPROXY_H_
#define RECRUITMENTLNODEPROXY_H_

class RecruitmentLnODEProxy  : public RecruitmentLnODE
{
private:
  virtual void solveODE(const int time, const Stats& statsPrevious, Stats& stats);

public:
  RecruitmentLnODEProxy();
  virtual ~RecruitmentLnODEProxy();
  /**
  * @copydoc GrSimulation::serialize
  */
  template<class Archive>
  void serialize(Archive& ar, const unsigned int version);

  RecruitmentMethod getMethod() const;

  RecruitmentBase* clone() const
  {
    return new RecruitmentLnODEProxy(*this);
  }

};

inline RecruitmentMethod RecruitmentLnODEProxy::getMethod() const
{
  return RECR_LN_ODE_PROXY;
}
template<class Archive>
inline void RecruitmentLnODEProxy::serialize(Archive& ar, const unsigned int /*version*/)
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(RecruitmentBase);
}

BOOST_CLASS_EXPORT_KEY(RecruitmentLnODEProxy)
#endif /* RECRUITMENTLNODEPROXY_H_ */
