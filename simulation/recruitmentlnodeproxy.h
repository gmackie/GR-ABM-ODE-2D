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

  RecruitmentMethod getMethod() const;
  void serialize(std::ostream&) const;
  void deserialize(std::istream&);

  RecruitmentBase* clone() const
  {
    return new RecruitmentLnODEProxy(*this);
  }

};

inline RecruitmentMethod RecruitmentLnODEProxy::getMethod() const
{
  return RECR_LN_ODE_PROXY;
}


// Even though this lymph ode proxy object uses the lymph node ode initial conditions array,
// it does not save any state between model timesteps. It stores values in the array based
// only on information for the current time step (ex. the time step value, time) and
// model parameters. So nothing needs to be done for serialization/deserialization.
inline void RecruitmentLnODEProxy::serialize(std::ostream&) const
{
}

inline void RecruitmentLnODEProxy::deserialize(std::istream&)
{
}

#endif /* RECRUITMENTLNODEPROXY_H_ */
