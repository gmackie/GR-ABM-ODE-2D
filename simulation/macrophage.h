/*
 * macrophage.h
 *
 *  Created on: 05-nov-2009
 *      Author: M. El-Kebir
 */

#ifndef MACROPHAGE_H
#define MACROPHAGE_H

#include "gr.h"
#include "agent.h"
#include "tgamma.h"

/*!
@brief Macrophage agent class.
@details
The state machine is given below:
\dot
digraph macrophage {
  edge [fontsize="10.0"];

  Resting [label="Resting, Mr" URL="\ref Mac::handleResting"];
  Infected [label="Infected, Mi" URL="\ref Mac::handleInfected"];
  Activated [label="Activated, Ma" URL="\ref Mac::handleActivated"];
  CInfected [label="Chronically Infected, Mci" URL="\ref Mac::handleCInfected"];
  DEAD [label="DEAD" URL="\ref Mac::handleDead"];

  Resting:n -> Resting:n [label="Killing of 1 Be IF (Be<=1)\nOR (Be>1 AND prob<25%)"];
  Resting -> Infected [label="Uptake of 1 Be IF Be>1 AND prob>25% then Bi=1"];
  Resting -> Activated [label="Mr activation induced by IFN-&gamma;\n(1 active T&gamma; cell in moore neighborhood)\n and TNF or Be"];
  Resting -> DEAD [label="Apoptosis (only TNF) OR Age"];

  Infected -> Infected [label="Uptake of Be IF (Be > 0 AND prob>(10-Bi)/100)"];
  Infected -> CInfected [label="Bi>10"];
  Infected -> Activated [label="Mi Activation (same as Mr)"];
  Infected -> DEAD [label="Apoptosis, Age OR CTL killing"];

  Activated -> DEAD [label="Apoptosis (only TNF) or Age"];

  CInfected -> DEAD [label="Apoptosis, Age, CTL killing OR bursting"];
}
\enddot
*/
class Mac : public Agent
{
public:
  enum State {MAC_RESTING, MAC_INFECTED, MAC_CINFECTED, MAC_ACTIVE, MAC_DEAD, NSTATES};
private:
  /*
   * !!! If the data members change then the serialize and deserialize functions need to be updated !!!
   */
  State _state;
  State _nextState;
  double _intMtb;
  bool _NFkB;
  bool _stat1;
  bool _ICOS;
  int _activationTime;
  int _deactivationTime;
  int _stat1Time;
  int _nfkbTime;

  /**
  * @brief
  * The probabilty for becoming infected is based on the following formula:
  * \f$\frac{Bi}{10} + \frac{prob}{1-25\%} = 1\f$
  * @param time
  * @param grid
  * @param stats
  * @param nfkbDynamics
  */
  void handleResting(const int time, GrGrid& grid, Stats& stats, bool nfkbDynamics);
  /**
  * @brief
  *
  * @param time
  * @param grid
  * @param stats
  * @param nfkbDynamics
  */
  void handleInfected(const int time, GrGrid& grid, Stats& stats, bool nfkbDynamics);
  /**
  * @brief
  *
  * @param time
  * @param grid
  * @param stats
  */
  void handleChronicallyInfected(const int time, GrGrid& grid, Stats& stats);
  /**
  * @details The following rules apply to recently activated macs:
  * -# STAT1 is activated (proxy for IFN-&gamma;) for (Mr)[Mac::handleResting] and (Mi)[Mac::handleInfected] with
  *    prob < 0.03 * (# T&gamma; in Moore neighborhood) - lifetime(=1.5 days)
  * -# AND NF-&kappa;B is activated: biologically this means either:
  *    + TNF > threshold_1 AND prob < (1-exp(boundTNFR1)
  *    + Bacteria (Be > 100, total in Moore neighborhood) - lifetime(=100minutes)
  * -# Ma kill internal bacteria (Bi) at a rate of 1 per time step (=10 minutes)
  *
  * @param time
  * @param grid
  * @param stats
  */
  void handleActivated(const int time, GrGrid& grid, Stats& stats);
  int getCountTgam(Tgam::State state, const GrGrid& grid) const;
  /**
  * @return Number of bacteria in moore neighborhood
  */
  double getExtMtbInMoore(const GrGrid& grid) const;

protected:
  friend class boost::serialization::access;
  Mac();

public:
  static double getIntMtbGrowthRate(const int time);

  Mac(int birthtime, int row, int col, Mac::State state, double intMtb, bool NFkB, bool stat1);
  ~Mac();

  static std::auto_ptr<ODESolvers::Stepper> stepper;
  static std::auto_ptr<LungFunc> deriv;

  /*virtual*/
  ODESolvers::Stepper* getStepper(ODESolvers::ODEMethod method)
  {
    static ODESolvers::ODEMethod initMethod = method;
    if(unlikely(Mac::stepper.get() == NULL || method != initMethod))
      {
        Mac::stepper.reset(ODESolvers::StepperFactory(method, _initvector.size()));
        initMethod = method;
      }
    return Mac::stepper.get();
  }

  static LungFunc* buildDerivFunc()
  {
    LungFunc* d = NULL;
    if(_PARAM(PARAM_NFKBODE_EN))
      d = new NFKBFunc(36);
    else
      d = Agent::buildDerivFunc();
    return d;
  }

  /*virtual*/ LungFunc* getDerivFunc()
  {
    if(deriv.get() == NULL)
      deriv.reset(Mac::buildDerivFunc());
    return deriv.get();
  }

  /*virtual*/ bool NFkBCapable() const
  {
    return true;
  }
  /*virtual*/ Agent* clone() const
  {
    return new Mac(*this);
  }
  void move(GrGrid& grid);
  /**
  * @copydoc Agent::secrete
  * The rules for secreting chemokines based on states and phenotype are as follows:
  * |     |              Nothing          |               NFkB              |            De-act             |
  * | --- | ----------------------------- | ------------------------------- | ----------------------------- |
  * | Mr  | No TNF<br/>No IL10<br/>No CCs | TNF 1/2<br/>No IL10<br/>CCs 1/2 | No TNF<br/>No IL10<br/>No CCs |
  * | Mi  | TNF 1/2<br/>IL10<br/>CCs 1/2  | TNF<br/>IL10<br/>CCs            | TNF 1/2<br/>IL10<br/>CCs 1/2  |
  * | Mci |                               | TNF<br/>2x IL10<br/>CCs         |                               |
  * | Ma  |                               | TNF<br/>IL10 1/10<br/>CCs       | No TNF<br/>No IL10<br/>No CCs |
  */
  void secrete(GrGrid& grid, bool tnfrDynamics, bool nfkbDynamics, bool tnfDepletion, bool il10rDynamics, bool il10Depletion, double mdt);
  void computeNextState(const int time, GrGrid& grid, Stats& stats, bool tnfrDynamics, bool nfkbDynamics, bool il10rDynamics, bool);
  void solveDegradation (GrGrid& grid, double dt, bool tnfrDynamics, bool il10rDynamics);
  void updateState();
  void updateStatistics(Stats& s) const;
  int getActivationTime() const;
  void setNFkB(bool value);
  bool getNFkB() const;
  bool getStat1() const;
  bool getICOS() const;
  int getState() const;
  Mac::State getNextState() const;
  double getIntMtb() const;
  void setIntMtb(double intMtb);
  void kill();
  void deactivate(const int time, Stats& stats);
  void apoptosis(GrGrid& grid);
  bool checkSTAT1(GrGrid& grid, const int time);
  bool checkNFkB(GrGrid& grid, const int time, bool tnfInducedNFkB);
  bool isDead() const;
  bool isDeadNext();
  static bool isMac(const Agent* pAgent);
  static bool isMac(const Agent* pAgent, Mac::State state);
  void print() const;
  template<class Archive>
  void serialize(Archive& ar, const unsigned int version);
  bool isDeactivated() const;
  double getNFkBn() const;
  double getNFkBc() const;
  double getNFkB_IkB() const;
  double getChemt() const;
  double getNormalizedACT() const;
  void setC1rrChemTNF(double value);
  AgentType getAgentType() const;
  void disperseMtb(GrGrid& grid, double fraction);
  template<typename Visitor>
  void visitProperties(Visitor& v);
  template<typename Visitor>
  void visitProperties(Visitor& v) const;
};

inline AgentType Mac::getAgentType() const
{
  return MAC;
}

inline bool Mac::isDeactivated() const
{
  return _deactivationTime != -1;
}

inline void Mac::setNFkB(bool value)
{
  _NFkB = value;
}

inline int Mac::getActivationTime() const
{
  return _activationTime;
}

inline bool Mac::getStat1() const
{
  return _stat1;
}

inline bool Mac::getNFkB() const
{
  return _NFkB;
}

inline bool Mac::getICOS() const
{
  return _ICOS;
}

inline int Mac::getState() const
{
  return (int)_state;
}

inline Mac::State Mac::getNextState() const
{
  return _nextState;
}

inline double Mac::getIntMtb() const
{
  return _intMtb;
}

inline void Mac::setIntMtb(double intMtb)
{
  _intMtb = intMtb;
}

inline bool Mac::isMac(const Agent* pAgent)
{
  return pAgent && pAgent->getAgentType() == MAC;
}

inline bool Mac::isMac(const Agent* pAgent, Mac::State state)
{
  if (!isMac(pAgent))
    {
      return false;
    }
  else
    {
      const Mac* pMac = static_cast<const Mac*>(pAgent);
      return pMac->getState() == state;
    }
}

inline bool Mac::isDead() const
{
  return _state == Mac::MAC_DEAD;
}

inline bool Mac::isDeadNext()
{
  return _nextState == Mac::MAC_DEAD;
}

inline double Mac::getNFkBn() const
{
  return _NFkBn;
}

inline double Mac::getNFkBc() const
{
  return _NFkBc;
}

inline double Mac::getNFkB_IkB() const
{
  return _NFkB_IkB;
}

inline double Mac::getChemt() const
{
  return _chemt;
}

inline double Mac::getNormalizedACT() const
{
  return _normalizedACT;
}

inline void Mac::setC1rrChemTNF(double value)
{
  _c1rrChemTNF = value;
}

inline double Mac::getIntMtbGrowthRate(const int time)
{
  double intMtbGrowthRate =  _PARAM(PARAM_INTMTB_GROWTH_RATE);

  if (time >= (_PARAM(PARAM_TCELL_TIME_RECRUITMENT_ENABLED) + _PARAM(PARAM_INTMTB_GROWTH_RATE_FACTOR_DELAY)))
    {
      // intMtb slows its growth rate as the immune system becomes activated. This is a proxy for more detailed intMtb growth rate response
      // to immune activation base on iNos and arginase.
      // The rate adjustment is applied to the fraction of the base growth rate, ex. for 1.025 it applies to the .025.
      // This ensures that the overall growth rate is positive. If applied to the entire base growth rate the effective
      // base growth rate can easily become < 1, i.e. negative growth.
      intMtbGrowthRate = 1 + (intMtbGrowthRate - 1) *_PARAM(PARAM_INTMTB_GROWTH_RATE_FACTOR_POST_ADAPTIVE);
    }

  return intMtbGrowthRate;
}

template<typename Visitor>
inline void Mac::visitProperties(Visitor& v)
{
  v.visit("intMtb", _intMtb, "");
  v.visit("NFkB", _NFkB, "");
  v.visit("stat1", _stat1, "");
  v.visit("ICOS", _ICOS, "");
  v.visit("activationTime", _activationTime, "");
  v.visit("deactivationTime", _deactivationTime, "");
  v.visit("stat1Time", _stat1Time, "");
  v.visit("nfkbTime", _nfkbTime, "");
  Agent::visitProperties(v);
}
template<typename Visitor>
inline void Mac::visitProperties(Visitor& v) const
{
  v.visit("intMtb", _intMtb, "");
  v.visit("NFkB", _NFkB, "");
  v.visit("stat1", _stat1, "");
  v.visit("ICOS", _ICOS, "");
  v.visit("activationTime", _activationTime, "");
  v.visit("deactivationTime", _deactivationTime, "");
  v.visit("stat1Time", _stat1Time, "");
  v.visit("nfkbTime", _nfkbTime, "");
  Agent::visitProperties(v);
}

std::ostream& operator<<(std::ostream& os, const Mac::State& s);


template<class Archive>
void Mac::serialize(Archive& ar, const unsigned int /*version*/) {
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Agent);
  ar & boost::serialization::make_nvp("state", _state);
  ar & boost::serialization::make_nvp("nextState", _nextState);
  ar & boost::serialization::make_nvp("intMtb", _intMtb);
  ar & boost::serialization::make_nvp("NFkB", _NFkB);
  ar & boost::serialization::make_nvp("stat1", _stat1);
  ar & boost::serialization::make_nvp("ICOS", _ICOS);
  ar & boost::serialization::make_nvp("activationTime", _activationTime);
  ar & boost::serialization::make_nvp("deactivationTime", _deactivationTime);
  ar & boost::serialization::make_nvp("stat1Time", _stat1Time);
  ar & boost::serialization::make_nvp("nfkbTime", _nfkbTime);
}

#endif /* MACROPHAGE_H */
