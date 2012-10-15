/*
 * recruitmentlnode.h
 *
 *  Created on: Jun 17, 2010
 *      Author: mohammed
 */

/*
 * This class originally was used to run a lymph node ode that used a DLL created by Matlab.
 * The Matlab DLL is no longer used, but everything in this function is used by subclasses,
 * except that they override the solveODE function.
 *
 */

#ifndef RECRUITMENTLNODE_H_
#define RECRUITMENTLNODE_H_

#include "gr.h"
#include "grgrid.h"
#include "recruitmentbase.h"
#include <string>

class RecruitmentLnODE : public RecruitmentBase
{
protected:

  static const int _nrConditions = 43;

  static const int _idxMDC = 13;
  static const int _idxNaiveCD4 = 14;
  static const int _idxNaiveCD8 = 16;
  static const int _idxEffectorTH1 = 39;
  static const int _idxEffectorT8 = 41;
  static const int _idxCTL = 42;

  double _tcellTable[TCELL_TYPE_COUNT]; // contains lower bounds

  int _tcellQueueCount[TCELL_TYPE_COUNT];
  typedef std::pair<int, TcellType> TcellTypePair;

  std::vector<TcellTypePair> _tcellQueue;

  int _prevMiMci;
  double _odeInitialConditions[_nrConditions];
  const std::string _odeApp;
  const std::string _odeTmpFile;

  void init();
  void updateInitialConditions(Stats& stats);
  virtual void solveODE(const int time, const Stats& statsPrevious, Stats& stats);
  void updateQueue(const int time, Stats& stats);
  void recruitMacsGetTcellSources(GrSimulation& sim, Stats& stats,
                                  ThresholdPosList tcellSources[TCELL_TYPE_COUNT]);
  void recruitMac(GrSimulation& sim, const Pos& pSource);
  void recruitTcells(GrSimulation& sim, Stats& stats,
                     ThresholdPosList tcellSources[TCELL_TYPE_COUNT]);

  friend class boost::serialization::access;
  RecruitmentLnODE() {}
public:
  RecruitmentLnODE(const std::string& odeApp, const std::string& odeTmpFile);
  virtual ~RecruitmentLnODE();

  RecruitmentMethod getMethod() const;
  template<class Archive>
  /**
  * @copydoc GrSimulation::serialize
  */
  void serialize(Archive& ar, const unsigned int version);
  virtual RecruitmentBase* clone() const
  {
    return new RecruitmentLnODE(*this);
  }

  void recruit(GrSimulation& sim, int time);
};

//@cond
// Register agent class as one that needs derived type translation
BOOST_SERIALIZATION_ASSUME_ABSTRACT(RecruitmentLnODE);
BOOST_CLASS_EXPORT_KEY(RecruitmentLnODE)
//@endcond

inline RecruitmentMethod RecruitmentLnODE::getMethod() const
{
  return RECR_LN_ODE;
}
template<class Archive>
inline void RecruitmentLnODE::serialize(Archive& ar, const unsigned int /*version*/)
{
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(RecruitmentBase);
  ar & boost::serialization::make_nvp("odeApp", const_cast<std::string&>(_odeApp));
  ar & boost::serialization::make_nvp("odeTmpFile", const_cast<std::string&>(_odeTmpFile));
  ar & BOOST_SERIALIZATION_NVP(_tcellTable);
  ar & BOOST_SERIALIZATION_NVP(_tcellQueueCount);
  ar & BOOST_SERIALIZATION_NVP(_tcellQueue);
  ar & BOOST_SERIALIZATION_NVP(_prevMiMci);
  ar & BOOST_SERIALIZATION_NVP(_odeInitialConditions);
}

#endif /* RECRUITMENTLNODE_H_ */
