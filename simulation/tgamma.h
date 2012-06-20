/*
 * tgamma.h
 *
 *  Created on: 03-nov-2009
 *      Author: M. El-Kebir
 */

#ifndef TGAMMA_H
#define TGAMMA_H

#include "gr.h"
#include "tcell.h"

class Tgam : public Tcell
{
public:
  enum State {TGAM_ACTIVE, TGAM_DOWN_REGULATED, TGAM_ACTIVE_DOUBLE, TGAM_INDUCED_REG, TGAM_DEAD, NSTATES};
private:
	static const std::string _ClassName;

	/*
	 * !!! If the data members change then the serialize and deserialize functions need to be updated !!!
	 */
	Tgam::State _state;
	Tgam::State _nextState;
	int _deactivationTime;
  int _transitionTime;
  
  int _nAntigenStim; // Number of successfull antigen stimulations
  int _nDownRegulated; // Number of downregulations
  int _nICOS; // Number of ICOS stimulations
	
	void handleActive(const int time, GrGrid& grid, Stats& stats, bool tgammatransition);
	void handleDownRegulated(const int time, GrGrid& grid, Stats& stats);
  void handleActiveDouble(const int time, GrGrid& grid, Stats& stats);
  void handleInducedReg(const int time, GrGrid& grid, Stats& stats);

public:
	Tgam();
	Tgam(int birthtime, int row, int col, Tgam::State state);
	~Tgam();
	void move(GrGrid& grid);
	void secrete(GrGrid& grid, bool tnfrDynamics, bool nfkbDynamics, bool tnfDepletion, bool il10rDynamics, bool il10Depletion, int mdt);
	void computeNextState(const int time, GrGrid& grid, Stats& stats, bool tnfrDynamics, bool nfkbDynamics, bool il10rDynamics, bool tgmmatransition);
	void updateState();
  void updateStatistics(Stats& s) const;
	int getState() const;
	Tgam::State getNextState() const;
	void deactivate(const int time);
	void kill();
	bool isDead() const;
	bool isDeadNext();
	int getDeactivationTime() const;
	static bool isTgam(const Agent* pAgent);
	static bool isTgam(const Agent* pAgent, Tgam::State state);
	void print() const;
	void serialize(std::ostream& out) const;
	void deserialize(std::istream& in);
	AgentType getAgentType() const;
  template<typename Visitor>
  void visitProperties(Visitor& v);
  template<typename Visitor>
  void visitProperties(Visitor& v) const;
};

inline AgentType Tgam::getAgentType() const
{
	return TGAM;
}

inline int Tgam::getDeactivationTime() const
{
	return _deactivationTime;
}

inline int Tgam::getState() const
{
	return _state;
}

inline Tgam::State Tgam::getNextState() const
{
	return _nextState;
}

inline bool Tgam::isTgam(const Agent* pAgent)
{
	return pAgent && pAgent->getAgentType() == TGAM;
}

inline bool Tgam::isTgam(const Agent* pAgent, Tgam::State state)
{
	if (!isTgam(pAgent))
	{
		return false;
	}
	else
	{
		const Tgam* pTgam = static_cast<const Tgam*>(pAgent);
		return pTgam->getState() == state;
	}
}

inline bool Tgam::isDead() const
{
	return _state == TGAM_DEAD;
}

inline bool Tgam::isDeadNext()
{
	return _nextState == TGAM_DEAD;
}

template<typename Visitor>
inline void Tgam::visitProperties(Visitor& v) {
  v.visit("deactivationTime", _deactivationTime, "");
  v.visit("transitionTime", _transitionTime, "");
  v.visit("AntigenStim", _nAntigenStim, "Number of successfull antigen stimulations");
  v.visit("DownRegulated", _nDownRegulated, "Number of downregulations");
  v.visit("ICOS", _nICOS, "Number of ICOS stimulations");
  Agent::visitProperties(v);
}
template<typename Visitor>
inline void Tgam::visitProperties(Visitor& v) const { 
  v.visit("deactivationTime", _deactivationTime, "");
  v.visit("transitionTime", _transitionTime, "");
  v.visit("AntigenStim", _nAntigenStim, "Number of successfull antigen stimulations");
  v.visit("DownRegulated", _nDownRegulated, "Number of downregulations");
  v.visit("ICOS", _nICOS, "Number of ICOS stimulations");
  Agent::visitProperties(v);
}

std::ostream& operator<<(std::ostream& os, const Tgam::State& s);

#endif /* TGAMMA_H */
