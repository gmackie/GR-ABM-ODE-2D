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

using namespace std;

class Mac : public Agent
{
public:
  enum State {MAC_RESTING, MAC_INFECTED, MAC_CINFECTED, MAC_ACTIVE, MAC_DEAD, NSTATES};
private:
	static const std::string _ClassName;
    static int _macodeSize;
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
	
	void handleResting(const int time, GrGrid& grid, Stats& stats, bool nfkbDynamics);
	void handleInfected(const int time, GrGrid& grid, Stats& stats, bool nfkbDynamics);
	void handleChronicallyInfected(const int time, GrGrid& grid, Stats& stats);
	void handleActivated(const int time, GrGrid& grid, Stats& stats);
	int getCountTgam(Tgam::State state, const GrGrid& grid) const;
	double getExtMtbInMoore(const GrGrid& grid) const;

public:
	static double getIntMtbGrowthRate(const int time);

	Mac();
	Mac(int birthtime, int row, int col, Mac::State state, double intMtb, bool NFkB, bool stat1);
	~Mac();

  static auto_ptr<ODESolvers::Stepper> stepper;
  static auto_ptr<ODESolvers::DerivativeFunc> deriv;

  /*virtual*/ ODESolvers::Stepper* getStepper(ODESolvers::ODEMethod method) {
    static ODESolvers::ODEMethod initMethod = method;
    if(unlikely(Mac::stepper.get() == NULL || method != initMethod)) {
      Mac::stepper.reset(ODESolvers::StepperFactory(method, _initvector.size()));
      initMethod = method;
    }
    return Mac::stepper.get();
  }

  static ODESolvers::DerivativeFunc* buildDerivFunc() {
    ODESolvers::DerivativeFunc* d = NULL;
    if(_PARAM(PARAM_NFKBODE_EN))
      d = new NFKBFunc(36);
    else
      d = Agent::buildDerivFunc();
    return d;
  }

  /*virtual*/ ODESolvers::DerivativeFunc* getDerivFunc() {
    if(deriv.get() == NULL)
      deriv.reset(Mac::buildDerivFunc());
    return deriv.get();
  }

  /*virtual*/ bool NFkBCapable() const { return true; }
  /*virtual*/ Agent* clone() const { return new Mac(*this); }
	void move(GrGrid& grid);
	void secrete(GrGrid& grid, bool tnfrDynamics, bool nfkbDynamics, bool tnfDepletion, bool il10rDynamics, bool il10Depletion, int mdt);
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
	void deactivate(const int time);
	void apoptosis(GrGrid& grid);
	bool isDead() const;
	bool isDeadNext();
	static bool isMac(const Agent* pAgent);
	static bool isMac(const Agent* pAgent, Mac::State state);
	void print() const;
	void serialize(std::ostream& out) const;
	void deserialize(std::istream& in);
	bool isDeactivated() const;
	double getNFkBn() const;
	double getNFkBc() const;
	double getNFkB_IkB() const;
	double getChemt() const;
	double getNormalizedACT() const;
	void setC1rrChemTNF(double value);
	AgentType getAgentType() const;
	void disperseMtb(GrGrid& grid, double fraction);
  static void setMacOdeSize(int odesize);
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

inline void Mac::setMacOdeSize(int odesize)
{
    _macodeSize = odesize;
}

template<typename Visitor>
inline void Mac::visitProperties(Visitor& v) {
  v.visit("intMtb", _intMtb, "");
  v.visit("NFkB", _NFkB, "");
  v.visit("stat1", _stat1, "");
  v.visit("ICOS", _ICOS, "");
  v.visit("activationTime", _activationTime, "");
  v.visit("deactivationTime", _deactivationTime, "");
  Agent::visitProperties(v);
}
template<typename Visitor>
inline void Mac::visitProperties(Visitor& v) const { 
  v.visit("intMtb", _intMtb, "");
  v.visit("NFkB", _NFkB, "");
  v.visit("stat1", _stat1, "");
  v.visit("ICOS", _ICOS, "");
  v.visit("activationTime", _activationTime, "");
  v.visit("deactivationTime", _deactivationTime, "");
  Agent::visitProperties(v);
}

std::ostream& operator<<(std::ostream& os, const Mac::State& s);

#endif /* MACROPHAGE_H */
