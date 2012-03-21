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

using namespace std;

class Mac : public Agent
{
private:
	static const std::string _ClassName;

	/*
	 * !!! If the data members change then the serialize and deserialize functions need to be updated !!!
	 */
	MacState _state;
	MacState _nextState;
	double _intMtb;
	bool _NFkB;
	bool _stat1;
	bool _ICOS;
	int _activationTime;
	int _deactivationTime;
    
	// NF-kB signaling pathway components
	double _IKKKa; // (IKKK in active state)
	double _IKKn; // (IKK in neutral state)
	double _IKKa; // (IKK in the active state)
	double _IKKi; // (IKK in inactive state)
	double _IkBp; // (Phospho-IkB)
	double _NFkB_IkBp; // %NFkB|IkBp
	double _NFkBc; // cytoplasmic NFkB
	double _NFkBn; // nucluar NFkB
	double _A20;
	double _A20t; // A20 transcript
	double _IkB;
	double _IkBn; // nucluar IkB
	double _IkBt; // IkB trancript
	double _NFkB_IkB;
	double _NFkB_IkBn;
	double _GA20; // (The state of A20)
	double _GIkB; // (The state of IkB)
	double _GR; // (The state of reporter genes)
	double _c1rrChemTNF; // NF-kB independent rate of TNF/chemokine mRNA synthesis  
	double _c1rChem;
	double _c1rTNF;
	double _chemt; // generic chemokine transcript
	double _chem; // intracellular generic chemokine protein
	double _TNFt; // TNF transcript
	double _TNF; // intracellular TNF 
	double _ACTt; // transcript of macrophage activating molecules
	double _ACT; // mac-activation molecules
	double _normalizedACT;
	double _IAPt; // transcript of IAP (inhibitor of apoptosis)
	double _IAP; 
	double _normalizedIAP;
	
	void handleResting(const int time, GrGrid& grid, GrStat& stats, bool nfkbDynamics);
	void handleInfected(const int time, GrGrid& grid, GrStat& stats, bool nfkbDynamics);
	void handleChronicallyInfected(const int time, GrGrid& grid, GrStat& stats);
	void handleActivated(const int time, GrGrid& grid, GrStat& stats);
	int getCountTgam(TgamState state, const GrGrid& grid) const;
	double getExtMtbInMoore(const GrGrid& grid) const;

public:
	static double getIntMtbGrowthRate(const int time);

	Mac();
	Mac(int birthtime, int row, int col, MacState state, double intMtb, bool NFkB, bool stat1);
	~Mac();
	void move(GrGrid& grid);
	void secrete(GrGrid& grid, bool tnfrDynamics, bool nfkbDynamics, bool tnfDepletion, bool il10rDynamics, bool il10Depletion);
	void computeNextState(const int time, GrGrid& grid, GrStat& stats, bool tnfrDynamics, bool nfkbDynamics, bool il10rDynamics, bool);
	void solveTNF (GrGrid& grid, double dt);
	void solveNFkBandTNF (GrGrid& grid, double dt);
	void solveNFkBODEsEquilibrium (double dt);
    void solveTNFandIL10 (GrGrid& grid, GrStat& stats, double dt, double currenttime);
    void solveTNFandIL10andNFkB (GrGrid& grid, double dt);
    void solveIL10 (GrGrid& grid, double dt);
    void solveDegradation (GrGrid& grid, double dt, bool tnfrDynamics, bool il10rDynamics);
	void updateState();
	int getActivationTime() const;
	void setNFkB(bool value);
	bool getNFkB() const;
	bool getStat1() const;
	bool getICOS() const;
	int getState() const;
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
	double getNFkBn() const;
	double getNFkBc() const;
	double getNFkB_IkB() const;
	double getChemt() const;
	double getNormalizedACT() const;
	void setC1rrChemTNF(double value);
	AgentType getAgentType() const;
	void disperseMtb(GrGrid& grid, double fraction);
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

inline bool Mac::getICOS() const
{
	return _ICOS;
}

inline int Mac::getState() const
{
	return (int)_state;
}

inline MacState Mac::getNextState() const
{
	return _nextState;
}

inline double Mac::getIntMtb() const
{
	return _intMtb;
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

inline double Mac::getNFkBn() const
{
	return _NFkBn;
}

inline double Mac::getNFkBc() const
{
	return _NFkBc;
}

inline double Mac::getNFkB_IkB() const
{
	return _NFkB_IkB;
}

inline double Mac::getChemt() const
{
	return _chemt;
}

inline double Mac::getNormalizedACT() const
{
	return _normalizedACT;
}

inline void Mac::setC1rrChemTNF(double value)
{
	_c1rrChemTNF = value;
}

inline double Mac::getIntMtbGrowthRate(const int time)
{
	double intMtbGrowthRate =  _PARAM(PARAM_INTMTB_GROWTH_RATE);

	if (time >= (_PARAM(PARAM_TCELL_TIME_RECRUITMENT_ENABLED) + _PARAM(PARAM_INTMTB_GROWTH_RATE_FACTOR_DELAY)))
	{
		// intMtb slows its growth rate as the immune system becomes activated. This is a proxy for more detailed intMtb growth rate response
		// to immune activation base on iNos and arginase.
		// The rate adjustment is applied to the fraction of the base growth rate, ex. for 1.025 it applies to the .025.
		// This ensures that the overall growth rate is positive. If applied to the entire base growth rate the effective
		// base growth rate can easily become < 1, i.e. negative growth.
		intMtbGrowthRate = 1 + (intMtbGrowthRate - 1) *_PARAM(PARAM_INTMTB_GROWTH_RATE_FACTOR_POST_ADAPTIVE);
	}

	return intMtbGrowthRate;
}

#endif /* MACROPHAGE_H */
