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

class Mac : public Agent
{
private:
	MacState _state;
	MacState _nextState;
	double _intMtb;
	bool _NFkB;
	bool _stat1;
	int _activationTime;
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
	
	void handleResting(const int time, GrGrid& grid, GrStat& stats);
	void handleInfected(const int time, GrGrid& grid, GrStat& stats);
	void handleChronicallyInfected(const int time, GrGrid& grid, GrStat& stats);
	void handleActivated(const int time, GrGrid& grid, GrStat& stats);
	int getCountTgam(TgamState state, const GrGrid& grid) const;
	double getExtMtbInMoore(const GrGrid& grid) const;

public:
	Mac(int birthtime, int row, int col, MacState state, double intMtb, bool NFkB, bool stat1);
	~Mac();
	void move(GrGrid& grid);
	void secrete(GrGrid& grid, bool tnfrDynamics);
	void computeNextState(const int time, GrGrid& grid, GrStat& stats, bool tnfrDynamics);
	void solveODEs (GrGrid& grid, double dt);
	void updateState();
	int getActivationTime() const;
	void setNFkB(bool value);
	bool getNFkB() const;
	bool getStat1() const;
	MacState getState() const;
	MacState getNextState() const;
	double getIntMtb() const;
	void setIntMtb(double intMtb);
	double getSurfBoundTNFR1() const;
	void kill();
	void deactivate(const int time);
	void apoptosis(GrGrid& grid);
	bool isDead();
	static bool isMac(const Agent* pAgent);
	static bool isMac(const Agent* pAgent, MacState state);
	void print() const;
	void serialize(std::ostream& out) const;
	void deserialize(std::istream& in);
	bool isDeactivated() const;
	AgentType getAgentType() const;
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

inline MacState Mac::getState() const
{
	return _state;
}

inline MacState Mac::getNextState() const
{
	return _nextState;
}

inline double Mac::getIntMtb() const
{
	return _intMtb;
}

inline double Mac::getSurfBoundTNFR1() const
{
	return _surfBoundTNFR1;
}

inline void Mac::setIntMtb(double intMtb)
{
	_intMtb = intMtb;
}

inline bool Mac::isMac(const Agent* pAgent)
{
	return pAgent && pAgent->getAgentType() == MAC;
}

inline bool Mac::isMac(const Agent* pAgent, MacState state)
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

inline bool Mac::isDead()
{
	return _state == MAC_DEAD;
}

#endif /* MACROPHAGE_H */
