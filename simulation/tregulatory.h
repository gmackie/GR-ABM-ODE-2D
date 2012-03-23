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
private:
	static const std::string _ClassName;

	/*
	 * !!! If the data members change then the serialize and deserialize functions need to be updated !!!
	 */
	TregState _state;
	TregState _nextState;
	void handleResting(const int time, GrGrid& grid, GrStat& stats);

public:
	Treg();
	Treg(int birthtime, int row, int col, TregState state);
	~Treg();
	void move(GrGrid& grid);
	void secrete(GrGrid& grid, bool tnfrDynamics, bool nfkbDynamics, bool tnfDepletion, bool il10rDynamics, bool il10Depletion);
	void computeNextState(const int time, GrGrid& grid, GrStat& stats, bool tnfrDynamics, bool nfkbDynamics, bool il10Depletion, bool);
	void updateState();
    void solveDegradation (GrGrid& grid, double dt, bool tnfrDynamics, bool il10rDynamics);
	int getState() const;
	TregState getNextState() const;
	void kill();
	void deactivate(const int time);
	bool isDead();
	static bool isTreg(const Agent* pAgent);
	static bool isTreg(const Agent* pAgent, TregState state);
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

inline TregState Treg::getNextState() const
{
	return _nextState;
}

inline bool Treg::isTreg(const Agent* pAgent)
{
	return pAgent && pAgent->getAgentType() == TREG;
}

inline bool Treg::isTreg(const Agent* pAgent, TregState state)
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

inline bool Treg::isDead()
{
	return _state == TREG_DEAD;
}

#endif /* TREGULATORY_H */
