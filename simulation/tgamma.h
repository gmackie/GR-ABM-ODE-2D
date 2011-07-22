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
	// TNF associated attributes
	double _mTNF; // No. of mTNF on the cell membrane
	double _surfTNFR1; // No. of cell surface TNFR1
	double _surfTNFR2;
	double _surfBoundTNFR1; // No. of sTNF-bound cell surface TNFR1
	double _surfBoundTNFR2;
	double _intBoundTNFR1; // No. of internalized TNF-bound TNFR1
	double _intBoundTNFR2;
	double _vTNFR1; // Rate of TNFR1 synthesis by cell
	double _vTNFR2;
	double _kSynth; // Rate of mTNF synthesis by cell
	double _kTACE; // Rate of mTNF release from cell by TACE activity
	
	void handleActive(const int time, GrGrid& grid, GrStat& stats);
	void handleDownRegulated(const int time, GrGrid& grid, GrStat& stats);

public:
	Tgam();
	Tgam(int birthtime, int row, int col, TgamState state);
	~Tgam();
	void move(GrGrid& grid);
	void secrete(GrGrid& grid, bool tnfrDynamics, bool nfkbDynamics, bool tnfDepletion);
	void computeNextState(const int time, GrGrid& grid, GrStat& stats, bool tnfrDynamics, bool nfkbDynamics);
	void updateState();
	void solveODEs (GrGrid& grid, double dt);
	TgamState getState() const;
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

inline TgamState Tgam::getState() const
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
