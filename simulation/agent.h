/*
 * gragent.h
 *
 *  Created on: 02-nov-2009
 *      Author: M. El-Kebir
 */

#ifndef AGENT_H
#define AGENT_H

#include "gr.h"
#include "grstat.h"
#include <string>

class Agent
{
protected:
	int _birthTime;
	int _deathTime;
	int _row;
	int _col;
	int moveAgent(GrGrid& grid, bool ccl2, bool ccl5, bool cxcl9, bool attractant, double bonusFactor);

public:
	Agent();
	Agent(int birthtime, int deathtime, int row, int col);
	virtual ~Agent();
	virtual void move(GrGrid& grid) = 0;
	virtual void secrete(GrGrid& grid, bool tnfrDynamics, bool tnfKnockout) = 0;
	virtual void computeNextState(const int time, GrGrid& grid, GrStat& stats, bool tnfrDynamics) = 0;
	virtual void updateState() = 0;
	virtual void kill() = 0;
	virtual void deactivate(const int time) = 0;
	virtual bool isDead() = 0;
	virtual void print() const = 0;
	virtual AgentType getAgentType() const = 0;
	int getRow() const;
	int getCol() const;
	int getBirthTime() const;
	int getDeathTime() const;
	bool timeToDie(const int time) const;
	virtual void serialize(std::ostream& out) const;
	virtual void deserialize(std::istream& in);
};

inline int Agent::getDeathTime() const
{
	return _deathTime;
}

inline bool Agent::timeToDie(const int time) const
{
	return time >= _deathTime;
}

inline int Agent::getBirthTime() const
{
	return _birthTime;
}

inline int Agent::getRow() const
{
	return _row;
}

inline int Agent::getCol() const
{
	return _col;
}

#endif /* AGENT_H */
