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
private:
	static const std::string _ClassName;

	/*
	 * !!! If the data members change then the serialize and deserialize functions need to be updated !!!
	 */
	TcytState _state;
	TcytState _nextState;
	int _deactivationTime;
	
	void handleActive(const int time, GrGrid& grid, GrStat& stats);
	void handleDownRegulated(const int time, GrGrid& grid, GrStat& stats);

public:
	Tcyt();
	Tcyt(int birthtime, int row, int col, TcytState state);
	~Tcyt();
	void move(GrGrid& grid);
	void secrete(GrGrid& grid, bool tnfrDynamics, bool nfkbDynamics, bool tnfDepletion, bool il10rDynamics, bool il10Depletion);
	void computeNextState(const int time, GrGrid& grid, GrStat& stats, bool tnfrDynamics, bool nfkbDynamics, bool il10rDynamics, bool);
	void updateState();
    void solveTNFandIL10 (GrGrid& grid, double dt);
    void solveIL10 (GrGrid& grid, double dt);
    void solveDegradation (GrGrid& grid, double dt, bool tnfrDynamics, bool il10rDynamics);
	int getState() const;
	TcytState getNextState() const;
	void deactivate(const int time);
	void kill();
	bool isDead();
	int getDeactivationTime() const;
	static bool isTcyt(const Agent* pAgent);
	static bool isTcyt(const Agent* pAgent, TcytState state);
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

inline TcytState Tcyt::getNextState() const
{
	return _nextState;
}

inline bool Tcyt::isTcyt(const Agent* pAgent)
{
	return pAgent && pAgent->getAgentType() == TCYT;
}

inline bool Tcyt::isTcyt(const Agent* pAgent, TcytState state)
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

#endif /* TCYTOTOXIC_H */
