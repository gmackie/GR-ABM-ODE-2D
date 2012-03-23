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
private:
	static const std::string _ClassName;

	/*
	 * !!! If the data members change then the serialize and deserialize functions need to be updated !!!
	 */
	TgamState _state;
	TgamState _nextState;
	int _deactivationTime;
    int _transitionTime;
    
    int _nAntigenStim; // Number of successfull antigen stimulations
    int _nDownRegulated; // Number of downregulations
    int _nICOS; // Number of ICOS stimulations
	
	void handleActive(const int time, GrGrid& grid, GrStat& stats, bool tgammatransition);
	void handleDownRegulated(const int time, GrGrid& grid, GrStat& stats);
    void handleActiveDouble(const int time, GrGrid& grid, GrStat& stats);
    void handleInducedReg(const int time, GrGrid& grid, GrStat& stats);

public:
	Tgam();
	Tgam(int birthtime, int row, int col, TgamState state);
	~Tgam();
	void move(GrGrid& grid);
	void secrete(GrGrid& grid, bool tnfrDynamics, bool nfkbDynamics, bool tnfDepletion, bool il10rDynamics, bool il10Depletion);
	void computeNextState(const int time, GrGrid& grid, GrStat& stats, bool tnfrDynamics, bool nfkbDynamics, bool il10rDynamics, bool tgmmatransition);
	void updateState();
	int getState() const;
	TgamState getNextState() const;
	void deactivate(const int time);
	void kill();
	bool isDead();
	int getDeactivationTime() const;
	static bool isTgam(const Agent* pAgent);
	static bool isTgam(const Agent* pAgent, TgamState state);
	void print() const;
	void serialize(std::ostream& out) const;
	void deserialize(std::istream& in);
	AgentType getAgentType() const;
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

inline TgamState Tgam::getNextState() const
{
	return _nextState;
}

inline bool Tgam::isTgam(const Agent* pAgent)
{
	return pAgent && pAgent->getAgentType() == TGAM;
}

inline bool Tgam::isTgam(const Agent* pAgent, TgamState state)
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

inline bool Tgam::isDead()
{
	return _state == TGAM_DEAD;
}

#endif /* TGAMMA_H */
