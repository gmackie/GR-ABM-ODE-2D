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
#include "numericalMethods.h"

using std::valarray;

struct LungFunc;

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
  
  #define AGENT_PROPS \
	P(unsigned long, id, 0, "Unique agent ID")  \
	P(int, birthTime, 0, "")  \
	P(int, deathTime, 0, "")  \
	P(Pos, pos, (Pos(-1, -1)), "")  \
	P(bool, trackMolecularDynamics, false, "")  \
	/* TNF associated attributes */ \
	P(Scalar, mTNF, 0.0, "No. of mTNF on the cell membrane")  \
	P(Scalar, surfTNFR1, 0.0, "No. of cell surface TNFR1") \
	P(Scalar, surfTNFR2, 0.0, "") \
	P(Scalar, surfBoundTNFR1, 0.0, "No. of sTNF-bound cell surface TNFR1") \
	P(Scalar, surfBoundTNFR2, 0.0, "") \
	P(Scalar, intBoundTNFR1, 0.0, "No. of internalized TNF-bound TNFR1") \
	P(Scalar, intBoundTNFR2, 0.0, "") \
  P(Scalar, mTNFRNA, 0.0, "") \
  P(Scalar, vTNFR1, 0.0, "Rate of TNFR1 synthesis by cell") \
  P(Scalar, vTNFR2, 0.0, "") \
  P(Scalar, kSynth, 0.0, "Rate of mTNF synthesis by cell") \
  P(Scalar, kTACE,  0.0, "Rate of mTNF release from cell by TACE activity") \
  P(Scalar, kmRNA,  0.0, "Rate of RNA synthesis for mTNF") \
  /* IL10 associated attributes */  \
  P(Scalar, surfIL10R, 0.0, "No. of cell surface IL10R") \
  P(Scalar, vIL10R, 0.0, "Rate of IL10R synthesis") \
  P(Scalar, surfBoundIL10R, 0.0, "No. of bound cell surface IL10R") \
  P(Scalar, kISynth, 0.0, "")  \
  P(Scalar, meanIL10R, 0.0, "") \
	/* NF-kB signaling pathway components  */ \
	P(Scalar, IKKKa, 0.0, "(IKKK in active state)") \
	P(Scalar, IKKn, (_PARAM(PARAM_GR_KNN)), "(IKK in neutral state)") \
	P(Scalar, IKKa, 0.0, "(IKK in the active state)") \
	P(Scalar, IKKi, 0.0, "(IKK in inactive state)") \
	P(Scalar, IkBp, 0.0, "(Phospho-IkB)") \
	P(Scalar, NFkB_IkBp, 0.0, "%NFkB|IkBp") \
	P(Scalar, NFkBc, 0.0, "cytoplasmic NFkB") \
	P(Scalar, NFkBn, 0.0, "nucluar NFkB") \
	P(Scalar, A20, 0.0, "") \
	P(Scalar, A20t, 0.0, "A20 transcript") \
	P(Scalar, IkB, 0.0, "") \
	P(Scalar, IkBn, 0.0, "nucluar IkB") \
	P(Scalar, IkBt, 0.0, "IkB trancript") \
	P(Scalar, NFkB_IkB, 0.0, "") \
	P(Scalar, NFkB_IkBn, 0.0, "") \
	P(Scalar, GA20, 0.0, "(The state of A20)") \
	P(Scalar, GIkB, 0.0, "(The state of IkB)") \
	P(Scalar, GR,   0.0, "(The state of reporter genes)") \
	P(Scalar, c1rrChemTNF, 0.0, "NF-kB independent rate of TNF/chemokine mRNA synthesis") \
	P(Scalar, c1rChem, 0.0, "") \
	P(Scalar, c1rTNF, 0.0, "") \
	P(Scalar, chemt, 0.0, "generic chemokine transcript") \
	P(Scalar, chem, 0.0, "intracellular generic chemokine protein") \
	P(Scalar, TNFt, 0.0, "TNF transcript") \
	P(Scalar, TNF, 0.0, "intracellular TNF") \
	P(Scalar, ACTt, 0.0, "transcript of macrophage activating molecules") \
	P(Scalar, ACT, 0.0, "mac-activation molecules") \
	P(Scalar, normalizedACT, 0.0, "") \
	P(Scalar, IAPt, 0.0, "transcript of IAP (inhibitor of apoptosis)") \
	P(Scalar, IAP, 0.0, "") \
	P(Scalar, normalizedIAP, 0.0, "")

#define P(type, name, ival, desc) type _##name; 
  AGENT_PROPS
#undef P

    bool _isODEsolved;
    Scalar _lasttimestep;
    std::valarray<double> _initvector;
    
	Pos moveAgent(GrGrid& grid, bool ccl2, bool ccl5, bool cxcl9, bool attractant, double bonusFactor);
	int getDestinationOrdinal(GrGrid& grid, bool ccl2, bool ccl5, bool cxcl9, bool attractant, double bonusFactor);
	Pos compartmentOrdinalToCoordinates(int ordinal, const Pos& dim) const;
  bool getODEstatus() const;
  void setODEstatus(bool isodesolved);

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

  static auto_ptr<ODESolvers::Stepper> stepper;
  static auto_ptr<ODESolvers::DerivativeFunc> deriv;

  virtual ODESolvers::Stepper* getStepper(ODESolvers::ODEMethod method) {
    static ODESolvers::ODEMethod initMethod = method;
    if(stepper.get() == NULL || method != initMethod) {
      stepper.reset(ODESolvers::StepperFactory(method, _initvector.size()));
      initMethod = method;
    }
    return stepper.get();
  }

  static ODESolvers::DerivativeFunc* buildDerivFunc(); 

  virtual ODESolvers::DerivativeFunc* getDerivFunc() {
    if(deriv.get() == NULL) {
      deriv.reset(buildDerivFunc());
    }
    return deriv.get();
  }

  virtual bool NFkBCapable() const { return false; }

	virtual void move(GrGrid& grid) = 0;
	virtual void secrete(GrGrid& grid, bool tnfrDynamics, bool nfkbDynamics, bool tnfDepletion, bool il10rDynamics, bool il10Depletion, int mdt) = 0;
	virtual void computeNextState(const int time, GrGrid& grid, Stats& stats, bool tnfrDynamics, bool nfkbDynamics, bool il10rDynamics, bool tgammatransition) = 0;
	virtual void updateState() = 0;
  virtual void updateStatistics(Stats& s) const = 0;
  void writeValarrayFromMembers(GrGrid& grid, valarray<double>& inputVector);
  void writeMembersFromValarray(GrGrid& grid, const valarray<double>& inputVector);

  void solveNFkBODEsEquilibrium (double dt);
  void solveMolecularScale(GrGrid& grid, double t, double dt, ODESolvers::ODEMethod method);
#if 0
    void solveTNF (GrGrid& grid, double dt);
    void solveIL10 (GrGrid& grid, double dt);
    void solveTNFandIL10(GrGrid& grid, double dt);
    void solveNFkBandTNF (GrGrid& grid, double dt);
    void solveTNFandIL10andNFkB (GrGrid& grid, double dt);
    
    void derivativeTNF(const valarray<double>& vecread, valarray<double>& vecwrite, double dt, GrGrid& grid);
    void derivativeIL10(const valarray<double>& vecread, valarray<double>& vecwrite, double dt, GrGrid& grid);
    void derivativeTNFandIL10(const valarray<double>& vecread, valarray<double>& vecwrite, double dt, GrGrid& grid);
    void derivativeTNFandNFKB(const valarray<double>& vecread, valarray<double>& vecwrite, double dt, GrGrid& grid);
    void derivativeTNFandIL10andNFKB(const valarray<double>& vecread, valarray<double>& vecwrite, double dt, GrGrid& grid);
    
    void solveRK2(GrGrid& grid, double dt, void(Agent::*derivativeType)(const valarray<double>& vecread, valarray<double>& vecwrite, double dt, GrGrid& grid));
    void solveRK4(GrGrid& grid, double dt, void(Agent::*derivativeType)(const valarray<double>& vecread, valarray<double>& vecwrite, double dt, GrGrid& grid));
    void solveForwardEuler(GrGrid& grid, double dt, void(Agent::*derivativeType)(const valarray<double>& vecread, valarray<double>& vecwrite, double dt, GrGrid& grid));
    void solveEulerPC(GrGrid& grid, double dt, void(Agent::*derivativeType)(const valarray<double>& vecread, valarray<double>& vecwrite, double dt, GrGrid& grid));

    void solveRKCK(GrGrid& grid, double dt, void(Agent::*derivativeType)(const valarray<double>& vecread, valarray<double>& vecwrite, double dt, GrGrid& grid));
    
    void RKStepper(GrGrid& grid, double dttry, double accuracy, double& dtnext, double& dtdid, double& currenttime, const double& starttime, const double& endtime, valarray<double> yscal, void(Agent::*derivativeType)(const valarray<double>& vecread, valarray<double>& vecwrite, double dt, GrGrid& grid));
    
    
    void AdaptiveRK(GrGrid& grid, double timestart, double timeend, double accuracy, double stepguess, double minstep, void(Agent::*derivs)(const valarray<double>& vecread, valarray<double>& vecwrite, double dt, GrGrid& grid));
    
    void RKStepper2cell(Agent* secondCell, GrGrid& grid, double dttry, double accuracy, double& dtnext, double& dtdid, double& currenttime, const double& starttime, const double& endtime, valarray<double> yscal, valarray<double> yscal2nd, void(Agent::*derivativeType)(const valarray<double>& vecread, valarray<double>& vecwrite, double dt, GrGrid& grid));
    
    
    void AdaptiveRK2cell(Agent* secondCell, GrGrid& grid, double timestart, double timeend, double accuracy, double stepguess, double minstep, void(Agent::*derivs)(const valarray<double>& vecread, valarray<double>& vecwrite, double dt, GrGrid& grid));
    
    void solveMolecularScaleFE(GrGrid& grid, double dt, bool nfkbDynamics, bool tnfrDynamics, bool il10rDynamics);
    void solveMolecularScaleRK2(GrGrid& grid, double dt, bool nfkbDynamics, bool tnfrDynamics, bool il10rDynamics);
    void solveMolecularScaleRK4(GrGrid& grid, double dt, bool nfkbDynamics, bool tnfrDynamics, bool il10rDynamics);
    void solveMolecularScaleRKadaptive(GrGrid& grid, double dt, bool nfkbDynamics, bool tnfrDynamics, bool il10rDynamics, double diffusionlength);
    void solveMolecularScaleEPC(GrGrid& grid, double dt, bool nfkbDynamics, bool tnfrDynamics, bool il10rDynamics);
#endif
  void checkTolerance(valarray<double>& veccheck);

	virtual void solveDegradation (GrGrid& grid, double dt, bool tnfrDynamics, bool il10rDynamics) = 0;
	virtual void solveDegradation (GrGrid& grid, double dt, bool tnfrDynamics, bool il10rDynamics, Scalar meanTNFR1, Scalar iIL10R);

	virtual void kill() = 0;
	virtual void deactivate(const int time) = 0;
	virtual bool isDead() const = 0;
	virtual bool isDeadNext() = 0;
	virtual void print() const = 0;
	virtual AgentType getAgentType() const = 0;
	virtual int getState() const = 0;
	int getRow() const;
	int getCol() const;
	bool timeToDie(const int time) const;
    
    virtual bool intCompareGT(const double param1, const double param2);
    

  #define P(type, name, ival, desc) \
    const type& get##name() const { return _##name; }
  AGENT_PROPS
  #undef P

	void setTrackMolecularDynamics(bool val) { _trackMolecularDynamics = val; }
  const Pos& getPosition() const { return _pos; } //Backwards compat
  unsigned long getID() const { return _id; }   //Backwards compat
	virtual void serialize(std::ostream& out) const;
	virtual void deserialize(std::istream& in);
  template<typename Visitor>
  void visitProperties(Visitor& v);
  template<typename Visitor>
  void visitProperties(Visitor& v) const;
};

inline unsigned long Agent::createID()
{
	return _nextID++;
}

inline bool Agent::timeToDie(const int time) const
{
	return time >= getdeathTime();
}

inline int Agent::getRow() const
{
	return _pos.x;
}

inline int Agent::getCol() const
{
	return _pos.y;
}

template<typename Visitor>
inline void Agent::visitProperties(Visitor& v) {
  #define P(type, name, ival, desc) \
    v.visit(#name, _##name, desc);
  AGENT_PROPS
  #undef P
}
template<typename Visitor>
inline void Agent::visitProperties(Visitor& v) const { 
  #define P(type, name, ival, desc) \
    v.visit(#name, _##name, desc);
  AGENT_PROPS
  #undef P
}
inline bool Agent::getODEstatus() const { return _isODEsolved; }
inline void Agent::setODEstatus(bool isodesolved) { _isODEsolved = isodesolved; }

struct LungFunc : ODESolvers::DerivativeFunc {
  const size_t il10offset;
  struct Params_t {
    Agent* agent;
    GrGrid* grid;
  };
  LungFunc(size_t il10) : DerivativeFunc(), il10offset(il10) {}
  virtual void operator()(const ODESolvers::ODEState& vecread, double /*t*/, ODESolvers::Derivative& vecwrite, void* params) const;
  void il10deriv(const ODESolvers::ODEState& vecread, double /*t*/, ODESolvers::Derivative& vecwrite, Params_t* params) const;
  virtual size_t dim() const { return _PARAM(PARAM_TNFODE_EN)*10+_PARAM(PARAM_IL10ODE_EN)*3; }
};

struct NFKBFunc : LungFunc {
  NFKBFunc(size_t il10) : LungFunc(il10) {}
  /*virtual*/ void operator()(const ODESolvers::ODEState& vecread, double /*t*/, ODESolvers::Derivative& vecwrite, void* params) const;
  /*virtual*/ size_t dim() const { return 36+_PARAM(PARAM_IL10ODE_EN)*3; }
};


inline ODESolvers::DerivativeFunc* Agent::buildDerivFunc() {
  LungFunc* d = new LungFunc(_PARAM(PARAM_TNFODE_EN)*10);
  return d;
}

#endif /* AGENT_H */
