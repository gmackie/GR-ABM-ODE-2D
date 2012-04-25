/*
 * tcytotoxic.h
 *
 *  Created on: 05-nov-2009
 *      Author: M. El-Kebir
 */

#ifndef TCYTOTOXIC_H
#define TCYTOTOXIC_H

#include "gr.h"
#include "tcell.h"

class Tcyt : public Tcell
{
public:
enum State {TCYT_DEAD, TCYT_ACTIVE, TCYT_DOWN_REGULATED, NSTATES};
private:
	static const std::string _ClassName;

	/*
	 * !!! If the data members change then the serialize and deserialize functions need to be updated !!!
	 */
	Tcyt::State _state;
	Tcyt::State _nextState;
	int _deactivationTime;
	
	void handleActive(const int time, GrGrid& grid, GrStat& stats);
	void handleDownRegulated(const int time, GrGrid& grid, GrStat& stats);

public:
	Tcyt();
	Tcyt(int birthtime, int row, int col, Tcyt::State state);
	~Tcyt();
	void move(GrGrid& grid);
	void secrete(GrGrid& grid, bool tnfrDynamics, bool nfkbDynamics, bool tnfDepletion, bool il10rDynamics, bool il10Depletion, int mdt);
	void computeNextState(const int time, GrGrid& grid, GrStat& stats, bool tnfrDynamics, bool nfkbDynamics, bool il10rDynamics, bool);
	void updateState();
	int getState() const;
	Tcyt::State getNextState() const;
	void deactivate(const int time);
	void kill();
	bool isDead();
	bool isDeadNext();
	int getDeactivationTime() const;
	static bool isTcyt(const Agent* pAgent);
	static bool isTcyt(const Agent* pAgent, Tcyt::State state);
	void print() const;
	void serialize(std::ostream& out) const;
	void deserialize(std::istream& in);
	AgentType getAgentType() const;
};

inline AgentType Tcyt::getAgentType() const
{
	return TCYT;
}

inline int Tcyt::getDeactivationTime() const
{
	return _deactivationTime;
}

inline int Tcyt::getState() const
{
	return _state;
}

inline Tcyt::State Tcyt::getNextState() const
{
	return _nextState;
}

inline bool Tcyt::isTcyt(const Agent* pAgent)
{
	return pAgent && pAgent->getAgentType() == TCYT;
}

inline bool Tcyt::isTcyt(const Agent* pAgent, Tcyt::State state)
{
	if (!isTcyt(pAgent))
	{
		return false;
	}
	else
	{
		const Tcyt* pTcyt = static_cast<const Tcyt*>(pAgent);
		return pTcyt->getState() == state;
	}
}

inline bool Tcyt::isDead()
{
	return _state == TCYT_DEAD;
}

inline bool Tcyt::isDeadNext()
{
	return _nextState == TCYT_DEAD;
}

#endif /* TCYTOTOXIC_H */
