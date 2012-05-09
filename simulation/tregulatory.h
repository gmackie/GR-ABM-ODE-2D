/*
 * tregulatory.h
 *
 *  Created on: 05-nov-2009
 *      Author: M. El-Kebir
 */

#ifndef TREGULATORY_H
#define TREGULATORY_H

#include "gr.h"
#include "tcell.h"

class Treg : public Tcell
{
public:
  enum State {TREG_ACTIVE, TREG_DEAD, NSTATES};
private:
	static const std::string _ClassName;

	/*
	 * !!! If the data members change then the serialize and deserialize functions need to be updated !!!
	 */
	Treg::State _state;
	Treg::State _nextState;
	void handleResting(const int time, GrGrid& grid, Stats& stats);

public:
	Treg();
	Treg(int birthtime, int row, int col, Treg::State state);
	~Treg();
	void move(GrGrid& grid);
	void secrete(GrGrid& grid, bool tnfrDynamics, bool nfkbDynamics, bool tnfDepletion, bool il10rDynamics, bool il10Depletion, int mdt);
	void computeNextState(const int time, GrGrid& grid, Stats& stats, bool tnfrDynamics, bool nfkbDynamics, bool il10Depletion, bool);
	void updateState();
  void updateStatistics(Stats& s) const;
	int getState() const;
	Treg::State getNextState() const;
	void kill();
	void deactivate(const int time);
	bool isDead() const;
	bool isDeadNext();
	static bool isTreg(const Agent* pAgent);
	static bool isTreg(const Agent* pAgent, Treg::State state);
	void print() const;
	void serialize(std::ostream& out) const;
	void deserialize(std::istream& in);
	AgentType getAgentType() const;
};

inline AgentType Treg::getAgentType() const
{
	return TREG;
}

inline int Treg::getState() const
{
	return _state;
}

inline Treg::State Treg::getNextState() const
{
	return _nextState;
}

inline bool Treg::isTreg(const Agent* pAgent)
{
	return pAgent && pAgent->getAgentType() == TREG;
}

inline bool Treg::isTreg(const Agent* pAgent, Treg::State state)
{
	if (!isTreg(pAgent))
	{
		return false;
	}
	else
	{
		const Treg* pTreg = static_cast<const Treg*>(pAgent);
		return pTreg->getState() == state;
	}
}

inline bool Treg::isDead() const
{
	return _state == TREG_DEAD;
}

inline bool Treg::isDeadNext()
{
	return _nextState == TREG_DEAD;
}

std::ostream& operator<<(std::ostream& os, const Treg::State& s);

#endif /* TREGULATORY_H */
