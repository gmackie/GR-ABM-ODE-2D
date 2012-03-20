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
private:
	static const std::string _ClassName;

	// A global ID counter for assigning a unique ID to an agent.
	static unsigned long nextID;

protected:
	/*
	 * !!! If the data members change then the serialize and deserialize functions need to be updated !!!
	 */
	static unsigned long createID();

	// Unique agent ID.
	unsigned long _id;
	int _birthTime;
	int _deathTime;
	Pos _pos;
	Pos moveAgent(GrGrid& grid, bool ccl2, bool ccl5, bool cxcl9, bool attractant, double bonusFactor);
	int getDestinationOrdinal(GrGrid& grid, bool ccl2, bool ccl5, bool cxcl9, bool attractant, double bonusFactor);
	Pos compartmentOrdinalToCoordinates(int ordinal, const Pos& dim) const;

public:
	Agent();
	Agent(int birthtime, int deathtime, int row, int col);
	virtual ~Agent();
	virtual void move(GrGrid& grid) = 0;
	virtual void secrete(GrGrid& grid, bool tnfrDynamics, bool nfkbDynamics, bool tnfDepletion, bool il10rDynamics, bool il10Depletion) = 0;
	virtual void computeNextState(const int time, GrGrid& grid, GrStat& stats, bool tnfrDynamics, bool nfkbDynamics, bool il10rDynamics, bool tgammatransition) = 0;
	virtual void updateState() = 0;
	virtual void kill() = 0;
	virtual void deactivate(const int time) = 0;
	virtual bool isDead() = 0;
	virtual void print() const = 0;
	virtual AgentType getAgentType() const = 0;
	virtual int getState() const = 0;
	unsigned long getID() const;
	const Pos& getPosition() const;
	int getRow() const;
	int getCol() const;
	int getBirthTime() const;
	int getDeathTime() const;
	bool timeToDie(const int time) const;
	virtual void serialize(std::ostream& out) const;
	virtual void deserialize(std::istream& in);
};

inline unsigned long Agent::createID()
{
	return nextID++;
}

inline unsigned long Agent::getID() const {
	return _id;
}

inline const Pos& Agent::getPosition() const {
  return _pos;
}

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
	return _pos.x;
}

inline int Agent::getCol() const
{
	return _pos.y;
}

#endif /* AGENT_H */
