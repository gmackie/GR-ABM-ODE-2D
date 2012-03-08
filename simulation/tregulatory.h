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
	// TNF associated attributes
	double _mTNF; // No. of mTNF on the cell membrane
	double _surfTNFR1; // No. of cell surface TNFR1
	double _surfTNFR2;
	double _surfBoundTNFR1; // No. of sTNF-bound cell surface TNFR1
	double _surfBoundTNFR2;
	double _intBoundTNFR1; // No. of internalized TNF-bound TNFR1
	double _intBoundTNFR2;
    double _mTNFRNA;
	double _vTNFR1; // Rate of TNFR1 synthesis by cell
	double _vTNFR2;
	double _kSynth; // Rate of mTNF synthesis by cell
	double _kTACE; // Rate of mTNF release from cell by TACE activity
    double _kmRNA; // Rate of RNA synthesis for mTNF
    
    // IL10 associated atributes
    double _surfIL10R; // No. of cell surface IL10R
    double _vIL10R; // Rate of IL10R synthesis
    double _surfBoundIL10R; // No. of bound cell surface IL10R
    double _kISynth;
	
	void handleResting(const int time, GrGrid& grid, GrStat& stats);

public:
	Treg();
	Treg(int birthtime, int row, int col, TregState state);
	~Treg();
	void move(GrGrid& grid);
	void secrete(GrGrid& grid, bool tnfrDynamics, bool nfkbDynamics, bool tnfDepletion, bool il10rDynamics, bool il10Depletion);
	void computeNextState(const int time, GrGrid& grid, GrStat& stats, bool tnfrDynamics, bool nfkbDynamics, bool il10Depletion, bool);
	void updateState();
	void solveTNF (GrGrid& grid, double dt);
    void solveTNFandIL10 (GrGrid& grid, double dt);
    void solveIL10 (GrGrid& grid, double dt);
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
