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
    double _kmRNA;
    
    // IL10 associated atributes
    double _surfIL10R; // No. of cell surface IL10R
    double _vIL10R; // Rate of IL10R synthesis
    double _surfBoundIL10R; // No. of bound cell surface IL10R
    double _kISynth;
    
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
	void solveTNF (GrGrid& grid, double dt);
    void solveTNFandIL10 (GrGrid&, double dt);
    void solveIL10 (GrGrid&, double dt);
    void solveDegradation (GrGrid& grid, double dt, bool tnfrDynamics, bool il10rDynamics);
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
