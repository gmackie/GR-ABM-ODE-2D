/*
 * gragent.h
 *
 *  Created on: 02-nov-2009
 *      Author: M. El-Kebir
 */

#ifndef AGENT_H
#define AGENT_H

#include "gr.h"
#include "params.h"
#include <string>
#include "numericalMethods.h"

using std::valarray;

struct LungFunc;

/**
* @brief Base agent class
* Abstract parent of all agent classes
* @note If the data members change then the serialize and deserialize functions need to be updated !!!
*/
class Agent
{
private:
  /// Fake RTTI for serialization purposes
	static const std::string _ClassName;

	/// A global ID counter for assigning a unique ID to an agent.
	static unsigned long _nextID;

protected:
  /**
  * @brief Creates an unique identifier for the agent
  *
  * @return a new id
  */
  static unsigned long createID();

  /**
  * @def AGENT_PROPS
  * X-Macro for holding all the agent-specific properties.
  * \see Agent::visitProperties Agent::serialize Agent::deserialize
  */
  
  /**
  * @def P(type, name, default, desc)
  * This macro is used to add a new property to an agent.  This will
  * automatically add code for serialization and accessing via corresponding
  * getname() and setname() functions.
  */
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
  \
  P(Scalar, M1M2Ratio, 0.0, "surfBoundTNFR1 / max(surfBoundIL10R1, 1) - updated in grsimulation") \
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

  /**
  * @name Agent Properties
  * Variable instantiations of all private properties of the agent class.
  * These variables are serialized and have accessor methods associated
  * @{
  */
#define P(type, name, ival, desc) type _##name;
  AGENT_PROPS
#undef P
  /**
  * @}
  */

  bool _isODEsolved;    /// Specifies whether this agent's odes have been solved for this timestep (for the two cell case)
  Scalar _lasttimestep; /// Last timestep used in an adaptive ode method
  std::valarray<double> _initvector;  /// valarray that holds all the integration variables to be passed to an integration method
    
  /**
  * @brief Move an agent based on chemotaxis
  *
  * @param grid Grid the agent is located on
  * @param ccl2 CCL2 effected?
  * @param ccl5 CCL5 effected?
  * @param cxcl9 CXCL9 effected?
  * @param attractant Mac Attractant effected?
  * @param bonusFactor Bonus weight for moving in the correct direction
  *
  * @return Calculated position to move to based on chemotaxis
  */
	Pos moveAgent(GrGrid& grid, bool ccl2, bool ccl5, bool cxcl9, bool attractant, double bonusFactor);
  /**
  * @brief returns an ordinal direction of the location to move based on chemotaxis
  *
  * @param grid Grid the agent is located on
  * @param ccl2 CCL2 effected?
  * @param ccl5 CCL5 effected?
  * @param cxcl9 CXCL9 effected?
  * @param attractant Mac Attractant effected?
  * @param bonusFactor Bonus weight for moving in the correct direction
  *
  * @return Calculated ordinial direction to move to based on chemotaxis
  * @see compartmentOrdinalToCoordinates
  */
	int getDestinationOrdinal(GrGrid& grid, bool ccl2, bool ccl5, bool cxcl9, bool attractant, double bonusFactor);
  /**
  * @brief Converts an ordinal to a grid position based on the agent's current position
  * The possible values for the ordinal are:
  *
  *  |   |   |   |
  *  | - | - | - |
  *  | 0 | 1 | 2 |
  *  | 3 | 4 | 5 |
  *  | 6 | 7 | 8 |
  *
  * @param ordinal Ordinal direction given by table above
  * @param dim Dimensions of the grid (to handle grid boundaries)
  *
  * @return Global position of indicated direction wrt to agent
  */
	Pos compartmentOrdinalToCoordinates(int ordinal, const Pos& dim) const;

public:
  /**
  * @brief Default constructor
  * Needed for deserializing the model state.
  * Avoids the calls to the random number generator in the normal constructor,
  * allowing the random number generator to remain in synch after
  * deserialization.
  */
	Agent();
  /**
  * @brief Standard constructor, must be used outside of serialization
  * @todo Finish documentation of these functions
  * @note This function uses the RNG to calculate internal variables.
  * @param birthtime
  * @param deathtime
  * @param row Initial row position of agent
  * @param col Initial column position of agent
  * @param meanTNFR1
  * @param stdTNFR1
  * @param meanTNFR2
  * @param stdTNFR2
  * @param kSynth
  * @param kTACE
  * @param iIL10R
  * @param stdIL10R
  * @param odesize number of ode variables (for setting up initial size of initvector)
  */
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

  /**
  * @brief serializes class specific meta-information
  * This function serializes meta-information (i.e. the last id used) to
  * maintain consistency between runs
  *
  * @param out output stream to write to
  */
	static void classSerialize(std::ostream& out);
  /**
  * @brief serializes class specific meta-information
  * This function serializes meta-information (i.e. the last id used) to
  * maintain consistency between runs
  *
  * @param in input stream to read from
  */
	static void classDeserialize(std::istream& in);

  static auto_ptr<ODESolvers::Stepper> stepper; /// ODE Stepper for all agents (unless overridden)
  static auto_ptr<LungFunc> deriv;              /// ODE Function for all agents (unless overridden)

  /**
  * @brief Gets the associated stepper for the agent type
  * Lazily builds the a stepper that can be used to integrate the ode function
  * with the associated method.
  *
  * @param method ODE method to use in order to integrate
  *
  * @return Stepper to use to integrate the method
  */
  virtual ODESolvers::Stepper* getStepper(ODESolvers::ODEMethod method) {
    static ODESolvers::ODEMethod initMethod = method;
    if(stepper.get() == NULL || method != initMethod) {
      stepper.reset(ODESolvers::StepperFactory(method, _initvector.size()));
      initMethod = method;
    }
    return stepper.get();
  }

  static LungFunc* buildDerivFunc(); 
  bool getODEstatus() const;
  void setODEstatus(bool isodesolved);

  virtual LungFunc* getDerivFunc() {
    if(deriv.get() == NULL) {
      deriv.reset(buildDerivFunc());
    }
    return deriv.get();
  }

  virtual bool NFkBCapable() const { return false; }

  /**
  * @brief Move the agent on the grid
  *
  * @param grid
  */
	virtual void move(GrGrid& grid) = 0;
  /**
  * @brief Causes the agent to secrete chemokines
  *
  * @param grid Grid that stores all chemokines
  * @param tnfrDynamics Account for tnfr dynamics?
  * @param nfkbDynamics Account for nfkb dynamics?
  * @param tnfDepletion Account for tnf Depletion?
  * @param il10rDynamics Account for il10r dynamics?
  * @param il10Depletion Account for il10 Depletion?
  * @param mdt Molecular timestep (how long to secrete for)
  */
  virtual void secrete(GrGrid& grid, bool tnfrDynamics, bool nfkbDynamics, bool tnfDepletion, bool il10rDynamics, bool il10Depletion, double mdt) = 0;
  /**
  * @brief Computes the next possible state for the agent
  * Given it's environment via the grid, computes the next state the agent should transition to for the next timestep
  * @param time Current time of the simulation (in seconds since infection)
  * @param grid Environment
  * @param stats Statistics to be updated
  * @param tnfrDynamics Account for tnfr dynamics
  * @param nfkbDynamics Account for nfkb dynamics
  * @param il10rDynamics Account for il10r dynamics
  * @param tgammatransition Account for T\gamma transitions
  */
	virtual void computeNextState(const int time, GrGrid& grid, Stats& stats, bool tnfrDynamics, bool nfkbDynamics, bool il10rDynamics, bool tgammatransition) = 0;
  /**
  * @brief Updates the current state of the agent to be the calculated next state \see computeNextState
  */
	virtual void updateState() = 0;
  /**
  * @brief Updates the statistics object with the final state of the agent
  *
  * @param s Statistics object to update
  */
  virtual void updateStatistics(Stats& s) const = 0;
  /**
  * @details Updates the M1M2Ratio based on it's internal SBTNFR1 and SBIL10R
  * using the following formula:
  * \f$M1M2ratio = \frac{IL10R1_{surfBound}}{\max (TNFR1_{surfBound}, 1)}\f$
  */
  void updateM1M2Ratio()
    { _M1M2Ratio = getsurfBoundTNFR1() / std::max(getsurfBoundIL10R(), 1.0); }
  /**
  * @brief Copies internal variables to valarray for solving ODEs
  *
  * @param grid
  * @param inputVector
  */
  void writeValarrayFromMembers(GrGrid& grid, valarray<double>& inputVector);
  /**
  * @brief Copies valarray ode variables to internal variables
  *
  * @param grid
  * @param inputVector
  */
  void writeMembersFromValarray(GrGrid& grid, const valarray<double>& inputVector);

  /**
  * @brief Solve internal ODEs for a preset time determined to give equilibrium
  *
  * @param dt Size of step to take in seconds
  */
  void solveNFkBODEsEquilibrium (double dt);
  /**
  * @brief Solves agent odes for 1 dt timestep using the given ode method
  *
  * @param grid
  * @param t Time of the simulation (in seconds since infection)
  * @param dt Time step size (in seconds)
  * @param method
  */
  void solveMolecularScale(GrGrid& grid, double t, double dt, ODESolvers::ODEMethod method);
  /**
  * @brief Solves agent odes for 1 dt timestep using the given ode method + adaptive algorithm
  * @note May be inaccurate due to coupled nature of multiple agents in a single grid space
  * @see adaptiveODE
  * @see adaptiveODE2Cell
  *
  * @param grid
  * @param t Time of the simulation (in seconds since infection)
  * @param dt Time step size (in seconds)
  * @param method
  */
  void solveMolecularScaleAdaptive(GrGrid& grid, double t, double dt, ODESolvers::ODEMethod method);
  void adaptiveODE(GrGrid& grid, double t, double dt, ODESolvers::ODEMethod method);
  void adaptiveODE2Cell(Agent* other, GrGrid& grid, double t, double dt, ODESolvers::ODEMethod method);

  bool TNFinducedApoptosis(GrGrid& grid, bool tnfrDynamics, bool nfkbDynamics);
  bool TNFinducedNFkB(GrGrid& grid, bool tnfrDynamics, bool nfkbDynamics);

  void checkTolerance(valarray<double>& veccheck);

	virtual void solveDegradation (GrGrid& grid, double dt, bool tnfrDynamics, bool il10rDynamics) = 0;
	virtual void solveDegradation (GrGrid& grid, double dt, bool tnfrDynamics, bool il10rDynamics, Scalar meanTNFR1, Scalar iIL10R);

  /**
  * @brief Kill the agent
  */
	virtual void kill() = 0;
  /**
  * @brief Deactivate the agent
  *
  * @param time Time (in seconds since infection) of deactivation
  * @param stats Statistics object to update
  */
  virtual void deactivate(const int time, Stats& stats) = 0;
  /**
  * @brief 
  *
  * @return True if dead, False otherwise
  */
	virtual bool isDead() const = 0;
  /**
  * @brief 
  *
  * @return True if the next state will be dead, False otherwise
  */
	virtual bool isDeadNext() = 0;
  /**
  * @brief Serializes the agent to stdout
  */
	virtual void print() const = 0;
  /**
  * @brief Returns the derived class's enumerated type tag for static casting
  *
  * @return Internal AgentType \see AgentType
  */
	virtual AgentType getAgentType() const = 0;
	virtual int getState() const = 0;
  /**
  * @brief Clones the agent
  * @return an exact duplicate of the agent
  */
  virtual Agent* clone() const = 0;
	int getRow() const;
	int getCol() const;
  /**
  * @brief 
  *
  * @param time Time (in seconds since infection) of simulation
  *
  * @return True if agent's life span has ended
  */
	bool timeToDie(const int time) const;
    
  /**
  * @brief Accurately compares two arguments to a specific precision
  *
  * @param param1
  * @param param2
  *
  * @return True if param1 > param2
  */
  virtual bool intCompareGT(const double param1, const double param2);
    

  /**
  * @name Property Accessors
  * @{
  */
  #define P(type, name, ival, desc) \
    const type& get##name() const { return _##name; }
  AGENT_PROPS
  #undef P
  /**@}*/

	void setTrackMolecularDynamics(bool val) { _trackMolecularDynamics = val; }
  const Pos& getPosition() const { return _pos; } //Backwards compat
  unsigned long getID() const { return _id; }   //Backwards compat
  /**
  * @brief Saves a serialized version of the agent to out \see deserialize
  *
  * @param out Output stream to save to
  */
	virtual void serialize(std::ostream& out) const;
  /**
  * @brief Loads an agent from the stream that was serialized \see serialize
  *
  * @param in Input stream to load from
  */
	virtual void deserialize(std::istream& in);
  /**
  * @brief visits all the properties of the agent class
  * @param v Visitor that will visit() each property
  * @tparam Visitor class that will visit each property via a visit(const char*, T&, const char*)
  */
  template<typename Visitor>
  void visitProperties(Visitor& v);
  /**
  * @brief visits all the properties of the agent class
  * @param v Visitor that will visit() each property
  * @tparam Visitor class that will visit each property via a visit(const char*, const T&, const char*)
  */
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
  virtual ~LungFunc() {}
  virtual void operator()(const ODESolvers::ODEState& vecread, double /*t*/, ODESolvers::Derivative& vecwrite, void* params) const;
  void il10deriv(const ODESolvers::ODEState& vecread, double /*t*/, ODESolvers::Derivative& vecwrite, Params_t* params) const;
  virtual size_t dim() const { return _PARAM(PARAM_TNFODE_EN)*10+_PARAM(PARAM_IL10ODE_EN)*3; }
  size_t tnfidx() const { return 8; }
  size_t il10idx() const { return il10offset; }
};

struct NFKBFunc : LungFunc {
  NFKBFunc(size_t il10) : LungFunc(il10) {}
  /*virtual*/ void operator()(const ODESolvers::ODEState& vecread, double /*t*/, ODESolvers::Derivative& vecwrite, void* params) const;
  /*virtual*/ size_t dim() const { return 36+_PARAM(PARAM_IL10ODE_EN)*3; }
};


inline LungFunc* Agent::buildDerivFunc() {
  LungFunc* d = new LungFunc(_PARAM(PARAM_TNFODE_EN)*10);
  return d;
}

#endif /* AGENT_H */
