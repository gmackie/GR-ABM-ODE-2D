/*
 * gragent.h
 *
 *  Created on: 02-nov-2009
 *      Author: M. El-Kebir
 */

#ifndef AGENT_H
#define AGENT_H

#include "gr.h"
#include <string>

using std::valarray;

class Agent
{
private:
	static const std::string _ClassName;

	// A global ID counter for assigning a unique ID to an agent.
	static unsigned long _nextID;

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

	bool _trackMolecularDynamics;

	// TNF associated attributes
	Scalar _mTNF; // No. of mTNF on the cell membrane
	Scalar _surfTNFR1; // No. of cell surface TNFR1
	Scalar _surfTNFR2;
	Scalar _surfBoundTNFR1; // No. of sTNF-bound cell surface TNFR1
	Scalar _surfBoundTNFR2;
	Scalar _intBoundTNFR1; // No. of internalized TNF-bound TNFR1
	Scalar _intBoundTNFR2;
    Scalar _mTNFRNA;
	Scalar _vTNFR1; // Rate of TNFR1 synthesis by cell
	Scalar _vTNFR2;
	Scalar _kSynth; // Rate of mTNF synthesis by cell
	Scalar _kTACE; // Rate of mTNF release from cell by TACE activity
    Scalar _kmRNA; // Rate of RNA synthesis for mTNF

    // IL10 associated attributes
    Scalar _surfIL10R; // No. of cell surface IL10R
    Scalar _vIL10R; // Rate of IL10R synthesis
    Scalar _surfBoundIL10R; // No. of bound cell surface IL10R
    Scalar _kISynth;
    Scalar _meanIL10R;
    
	// NF-kB signaling pathway components
	Scalar _IKKKa; // (IKKK in active state)
	Scalar _IKKn; // (IKK in neutral state)
	Scalar _IKKa; // (IKK in the active state)
	Scalar _IKKi; // (IKK in inactive state)
	Scalar _IkBp; // (Phospho-IkB)
	Scalar _NFkB_IkBp; // %NFkB|IkBp
	Scalar _NFkBc; // cytoplasmic NFkB
	Scalar _NFkBn; // nucluar NFkB
	Scalar _A20;
	Scalar _A20t; // A20 transcript
	Scalar _IkB;
	Scalar _IkBn; // nucluar IkB
	Scalar _IkBt; // IkB trancript
	Scalar _NFkB_IkB;
	Scalar _NFkB_IkBn;
	Scalar _GA20; // (The state of A20)
	Scalar _GIkB; // (The state of IkB)
	Scalar _GR; // (The state of reporter genes)
	Scalar _c1rrChemTNF; // NF-kB independent rate of TNF/chemokine mRNA synthesis  
	Scalar _c1rChem;
	Scalar _c1rTNF;
	Scalar _chemt; // generic chemokine transcript
	Scalar _chem; // intracellular generic chemokine protein
	Scalar _TNFt; // TNF transcript
	Scalar _TNF; // intracellular TNF 
	Scalar _ACTt; // transcript of macrophage activating molecules
	Scalar _ACT; // mac-activation molecules
	Scalar _normalizedACT;
	Scalar _IAPt; // transcript of IAP (inhibitor of apoptosis)
	Scalar _IAP; 
	Scalar _normalizedIAP;

    valarray<double> _initvector;
    valarray<double> _k1vector;
    valarray<double> _k2vector;
    valarray<double> _k3vector;
    valarray<double> _k4vector;
    valarray<double> _switchvector;
    
	Pos moveAgent(GrGrid& grid, bool ccl2, bool ccl5, bool cxcl9, bool attractant, double bonusFactor);
	int getDestinationOrdinal(GrGrid& grid, bool ccl2, bool ccl5, bool cxcl9, bool attractant, double bonusFactor);
	Pos compartmentOrdinalToCoordinates(int ordinal, const Pos& dim) const;

public:
	Agent();
	Agent(int birthtime, int deathtime, int row, int col
			//TNFR Components
			, Scalar meanTNFR1
			, Scalar stdTNFR1
			, Scalar meanTNFR2
			, Scalar stdTNFR2
			, Scalar kSynth
			, Scalar kTACE

			// IL10 components
			, Scalar iIL10R
			, Scalar stdIL10R
            , int odesize
		);

	virtual ~Agent();

	static void classSerialize(std::ostream& out);
	static void classDeserialize(std::istream& in);

	virtual void move(GrGrid& grid) = 0;
	virtual void secrete(GrGrid& grid, bool tnfrDynamics, bool nfkbDynamics, bool tnfDepletion, bool il10rDynamics, bool il10Depletion, int mdt) = 0;
	virtual void computeNextState(const int time, GrGrid& grid, GrStat& stats, bool tnfrDynamics, bool nfkbDynamics, bool il10rDynamics, bool tgammatransition) = 0;
	virtual void updateState() = 0;

	virtual void solveTNF (GrGrid& grid, double dt);
	virtual void solveIL10 (GrGrid& grid, double dt);
	virtual void solveTNFandIL10(GrGrid& grid, double dt);
    virtual void solveNFkBandTNF (GrGrid& grid, double dt);
    virtual void solveTNFandIL10andNFkB (GrGrid& grid, double dt);
    virtual void solveNFkBODEsEquilibrium (double dt);
    
    virtual void derivativeTNF(const valarray<double>& vecread, valarray<double>& vecwrite, double dt, GrGrid& grid);
    virtual void derivativeIL10(const valarray<double>& vecread, valarray<double>& vecwrite, double dt, GrGrid& grid);
    virtual void derivativeTNFandIL10(const valarray<double>& vecread, valarray<double>& vecwrite, double dt, GrGrid& grid);
    virtual void derivativeTNFandNFKB(const valarray<double>& vecread, valarray<double>& vecwrite, double dt, GrGrid& grid);
    virtual void derivativeTNFandIL10andNFKB(const valarray<double>& vecread, valarray<double>& vecwrite, double dt, GrGrid& grid);
    virtual void writeValarrayFromMembers(GrGrid& grid, valarray<double>& inputVector);
    virtual void writeMembersFromValarray(GrGrid& grid, const valarray<double>& inputVector);
    
    virtual void solveRK2(GrGrid& grid, double dt, void(Agent::*derivativeType)(const valarray<double>& vecread, valarray<double>& vecwrite, double dt, GrGrid& grid));
    virtual void solveRK4(GrGrid& grid, double dt, void(Agent::*derivativeType)(const valarray<double>& vecread, valarray<double>& vecwrite, double dt, GrGrid& grid));
    virtual void solveForwardEuler(GrGrid& grid, double dt, void(Agent::*derivativeType)(const valarray<double>& vecread, valarray<double>& vecwrite, double dt, GrGrid& grid));
    virtual void solveEulerPC(GrGrid& grid, double dt, void(Agent::*derivativeType)(const valarray<double>& vecread, valarray<double>& vecwrite, double dt, GrGrid& grid));

    
    virtual void solveMolecularScaleFE(GrGrid& grid, double dt, bool nfkbDynamics, bool tnfrDynamics, bool il10rDynamics);
    virtual void solveMolecularScaleRK2(GrGrid& grid, double dt, bool nfkbDynamics, bool tnfrDynamics, bool il10rDynamics);
    virtual void solveMolecularScaleRK4(GrGrid& grid, double dt, bool nfkbDynamics, bool tnfrDynamics, bool il10rDynamics);
    virtual void solveMolecularScaleEPC(GrGrid& grid, double dt, bool nfkbDynamics, bool tnfrDynamics, bool il10rDynamics);


    
	virtual void solveDegradation (GrGrid& grid, double dt, bool tnfrDynamics, bool il10rDynamics) = 0;
	virtual void solveDegradation (GrGrid& grid, double dt, bool tnfrDynamics, bool il10rDynamics, Scalar meanTNFR1, Scalar iIL10R);

	virtual void kill() = 0;
	virtual void deactivate(const int time) = 0;
	virtual bool isDead() = 0;
	virtual bool isDeadNext() = 0;
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

	// TNF associated attributes
	Scalar getMTNF() const;
	Scalar getSurfTNFR1() const;
	Scalar getSurfTNFR2() const;
	Scalar getSurfBoundTNFR1() const;
	Scalar getSurfBoundTNFR2() const;
	Scalar getIntBoundTNFR1() const;
	Scalar getIntBoundTNFR2() const;
	Scalar getMTNFRNA() const;
	Scalar getVTNFR1() const;
	Scalar getVTNFR2() const;
	Scalar getKSynth() const;
	Scalar getKTACE() const;
	Scalar getKmRNA() const;

	// IL10 associated attributes
	Scalar getSurfIL10R() const;
	Scalar getVIL10R() const;
	Scalar getSurfBoundIL10R() const;
	Scalar getKISynth() const;


	bool getTrackMolecularDynamics() const;
	void setTrackMolecularDynamics(bool trackMolecularDynamics);
	virtual void serialize(std::ostream& out) const;
	virtual void deserialize(std::istream& in);
};

inline unsigned long Agent::createID()
{
	return _nextID++;
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

// TNF associated attributes
inline Scalar Agent::getMTNF() const
{
	 return _mTNF;
}

inline Scalar Agent::getSurfTNFR1() const
{
	 return _surfTNFR1;
}

inline Scalar Agent::getSurfTNFR2() const
{
	 return _surfTNFR2;
}

inline Scalar Agent::getSurfBoundTNFR1() const
{
	return _surfBoundTNFR1;
}

inline Scalar Agent::getSurfBoundTNFR2() const
{
	 return _surfBoundTNFR2;
}

inline Scalar Agent::getIntBoundTNFR1() const
{
	 return _intBoundTNFR1;
}

inline Scalar Agent::getIntBoundTNFR2() const
{
	 return _intBoundTNFR2;
}

inline Scalar Agent::getMTNFRNA() const
{
	 return _mTNFRNA;
}

inline Scalar Agent::getVTNFR1() const
{
	 return _vTNFR1;
}

inline Scalar Agent::getVTNFR2() const
{
	 return _vTNFR2;
}

inline Scalar Agent::getKSynth() const
{
	 return _kSynth;
}

inline Scalar Agent::getKTACE() const
{
	 return _kTACE;
}

inline Scalar Agent::getKmRNA() const
{
	 return _kmRNA;
}

// IL10 associated attributes
inline Scalar Agent::getSurfIL10R() const
{
	 return _surfIL10R;
}

inline Scalar Agent::getVIL10R() const
{
	 return _vIL10R;
}

inline Scalar Agent::getSurfBoundIL10R() const
{
	 return _surfBoundIL10R;
}

inline Scalar Agent::getKISynth() const
{
	 return _kISynth;
}

inline bool Agent::getTrackMolecularDynamics() const
{
	return _trackMolecularDynamics;
}

inline void Agent::setTrackMolecularDynamics(bool trackMolecularDynamics)
{
	_trackMolecularDynamics = trackMolecularDynamics;
}

#endif /* AGENT_H */
