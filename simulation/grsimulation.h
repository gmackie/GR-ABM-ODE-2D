/*
 * grsimulation.h
 *
 *  Created on: 05-nov-2009
 *      Author: M. El-Kebirƒ
 */

#ifndef GRSIMULATION_H
#define GRSIMULATION_H

#include "gr.h"
#include "grsimulationgrid.h"
#include "stat.h"
#include "macrophage.h"
#include "tgamma.h"
#include "tcytotoxic.h"
#include "tregulatory.h"
#include "ttest.h"
#include "grdiffusion.h"
#include "recruitmentbase.h"
#include "params.h"
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/placeholders.hpp>

/**
* @brief Simulation class that holds all information about the simulation.
* This includes interacting agents, simulation options, diffusing chemicals,
* etc.
*
* @note If data members change, serialization must also change to be consistent
*/
class GrSimulation
{
private:

  /// @cond
  /// Functors for boost serialization
  template<typename T> struct wrap {};
  template<typename Archive> struct serialization_class_register_t {
    Archive& ar;
    serialization_class_register_t(Archive& _ar) : ar(_ar) {}
    template<typename T>
    void operator()(wrap<T>)
      { ar.template register_type<T>(); }
  };
  /// @endcond

  int _time;
  GrSimulationGrid _grid;
  /**
  * @name Agent Lists
  * @TODO Replace this with a templated(?) map in order to decouple agent type from grsimulation
  * @brief These lists hold pointers to each agent type.
  * @{
  */
  MacList _macList;
  TgamList _tgamList;
  TcytList _tcytList;
  TregList _tregList;
  /**@}*/

  Stats _statsPrevious; /// The stats at the start of a time step - from the end of the previous time step.
  Stats _stats;

  /// @TODO Documentation
  double _areaThreshold;
  /// @TODO Documentation
  double _areaThresholdCellDensity;
  /// Diffusion method to use
  GrDiffusion* _pDiffusion;
  /// @TODO Documentation
  TTest* _pTTest[NOUTCOMES];
  /// Recruitment method to use
  RecruitmentBase* _pRecruitment;
  /// Enable TNFR Dynamics?
  bool _tnfrDynamics;
  /// Enable NFkB Dynamics?
  bool _nfkbDynamics;
  /// Enable IL10R Dynamics?
  bool _il10rDynamics;
  /// @TODO Documentation
  bool _tgammatransition;
  /// Use the adaptive method?
  bool _adaptive;

  /// ODE Solver method to use for agent odes
  ODESolvers::ODEMethod _odeSolver;

  /// Inhibits tnf secretion if true and if not using tnfr dynamics.
  int _tnfDepletionTimeStep;
  /// Inhibits il10 secretion if true and if not using il10r dynamics.
  int _il10DepletionTimeStep;

  /**
  * Whether or not TCell recruitment has been enabled.
  * Once enabled it stays enabled even if the criteria by which it became enabled changes.
  */
  bool _tcellRecruitmentBegun;
  /// Number of timesteps of mdt (see \ref Params) for diffusion
  int _numMolecularPerDiffusion;
  /// Number of timesteps of dt (see \ref Params) for agent odes
  int _numDiffusionPerAgent;

  void moveTcells();
  void moveMacrophages();
  void updateStates();
  void updateT_Test();
  void computeNextStates();
  void secreteFromMacrophages(bool tnfDepletion, bool il10Depletion, double mdt);
  void secreteFromTcells(bool tnfDepletion, bool il10Depletion, double mdt);
  void secreteFromCaseations(double mdt);
  void solveMolecularScale(double dt);
  void solveMolecularScaleAdaptive(double dt);
  void adjustTNFDegradation(double dt);
  void adjustFauxDegradation(double dt);
  void growExtMtb();
  void shuffleCells();
  void checkTCellRecruitmentStart();

public:
  /**
  * @param dim Dimensions of the internal grid to use (in # of 20 micron sectons)
  */
  GrSimulation(const Pos& dim);
  ~GrSimulation();
  /**
  * @brief Initialize the simulation for t=0 with initial macs, etc
  */
  void init();
  /**
  * @brief Select all macrophages within the radius to be used in tracking
  *
  * @param molecularTrackingRadius defined radius (must be 0<r<DIM)
  */
  void initMolecularTracking(Scalar molecularTrackingRadius);
  /**
  * @brief Select all macrophages given to be used in tracking
  *
  * @param ids list of agent ids
  */
  void initMolecularTracking(const std::vector<size_t>& ids);
  /**
  * @brief Run the simulation for a 10 minute timestep
  */
  void solve();
  void performT_Test();
  /**
  * @name Accessors
  * @{
  */

  /**
  * @return Number of 10 minute timesteps since infection
  */
  int getTime() const;

  const Stats& getStats() const;
  Stats& getStats();

  const Stats& getStatsPrevious() const;
  Stats& getStatsPrevious();

  const GrGrid& getGrid() const;
  GrGrid& getGrid();
  const MacList& getMacList() const;
  const TgamList& getTgamList() const;
  const TcytList& getTcytList() const;
  const TregList& getTregList() const;
  DiffusionMethod getDiffusionMethod() const;
  void setDiffusionMethod(DiffusionMethod method);

  RecruitmentMethod getRecruitmentMethod();
  void setRecruitmentMethod(RecruitmentMethod method);

  void setODESolverMethod(ODESolvers::ODEMethod odemethod);
  int getODESolverMethod() const;
  bool getTnfrDynamics() const;
  void setTnfrDynamics(bool tnfrDynamics);
  bool getNfkbDynamics() const;
  void setNfkbDynamics(bool nfkbDynamics);
  int getTnfDepletionTimeStep() const;
  void setTnfDepletionTimeStep(int tnfDepletionTimeStep);
  bool getIl10rDynamics() const;
  void setIl10rDynamics(bool il10rdynamics);
  int getIl10DepletionTimeStep() const;
  void setIl10DepletionTimeStep (int il10DepletionTimeStep);
  bool getTgammaTransition() const;
  void setTgammaTransition(bool tgammatransition);
  bool getAdaptive() const;
  void setAdaptive(bool adaptive);
  double getAreaThreshold() const;
  void setAreaThreshold(double areaThreshold);
  double getAreaThresholdCellDensity() const;
  void setAreaThresholdCellDensity(double areaThreshold);
  OutcomeMethod getOutcomeMethod(int index) const;
  void getOutcomeParameters(int index, int& samplePeriod, int& testPeriod, double& alpha) const;
  void setOutcomeMethod(int index, OutcomeMethod method, double alpha, int testPeriod, int samplePeriod);
  bool getTCellRecruitmentBegun();
  /**@}*/

  /**
  * Convenience function for saving the simulation state
  * @see load
  * @see serialize
  */
  void save(const char* fname) const;
  /// @overload void GrSimulation::save(const char* fname)
  void save(std::ostream& out) const;
  /**
  * @brief 
  *
  * @tparam Archive 
  * @param ar boost::archive to save the object to.
  * @param version version number of the class to be serialized (if SVN_VERSION defined, this will be the revision number)
  */
  template<class Archive>
  void serialize(Archive& ar, const unsigned int version);
  /**
  * Convenience function for loading the simulation state
  * @see save
  * @see serialize
  */
  void load(const char* fname);
  /// @overload void GrSimulation::load(const char* fname)
  void load(std::istream& in);
  void deserializeRecruitmentMethod(RecruitmentMethod method, std::istream& in);
  /**
  * @name Agent Factory Methods
  * @todo Replace with a templated(?) form to decouple agent type
  * @{
  */
  Mac* createMac(int row, int col, int birthtime, Mac::State state, bool NFkB, bool stat1);
  Tgam* createTgam(int row, int col, int birthtime, Tgam::State state);
  Tcyt* createTcyt(int row, int col, int birthtime, Tcyt::State state);
  Treg* createTreg(int row, int col, int birthtime, Treg::State state);
  /**@}*/

  /**
  * @brief Converts the time given by getTime() to appropriate days, hours, minutes
  *
  * @param[in] time
  * @param[out] rDays
  * @param[out] rHours
  * @param[out] rMinutes
  */
  static void convertSimTime(const int time, int& rDays, int& rHours, int& rMinutes);
  void timestepSync();
  /**
  * @brief Deep-copy of the simulation
  */
  GrSimulation* clone() const;
};

inline bool GrSimulation::getTCellRecruitmentBegun()
{
  return _tcellRecruitmentBegun;
}

inline void GrSimulation::getOutcomeParameters(int index,
    int& samplePeriod, int& testPeriod, double& alpha) const
{
  assert(0 <= index && index < NOUTCOMES);

  if (_pTTest[index])
    {
      samplePeriod = _pTTest[index]->getSamplePeriod();
      testPeriod = _pTTest[index]->getTestPeriod();
      alpha = _pTTest[index]->getAlpha();
    }
  else
    {
      samplePeriod = 1;
      testPeriod = 2;
      alpha = 0.05;
    }
}

inline OutcomeMethod GrSimulation::getOutcomeMethod(int index) const
{
  assert(0 <= index && index < NOUTCOMES);
  return _pTTest[index] ? _pTTest[index]->getMethod() : OUTCOME_NONE;
}

inline double GrSimulation::getAreaThreshold() const
{
  return _areaThreshold;
}

inline void GrSimulation::setAreaThreshold(double areaThreshold)
{
  _areaThreshold = areaThreshold;
}

inline double GrSimulation::getAreaThresholdCellDensity() const
{
  return _areaThresholdCellDensity;
}

inline void GrSimulation::setAreaThresholdCellDensity(double areaThreshold)
{
  _areaThresholdCellDensity = areaThreshold;
}

inline void GrSimulation::convertSimTime(const int time, int& rDays, int& rHours, int& rMinutes)
{
  rDays = time / TIME_STEPS_PER_DAY;
  rHours = (time % TIME_STEPS_PER_DAY) / 6;
  rMinutes = ((time % TIME_STEPS_PER_DAY) % 6) * 10;
}

inline int GrSimulation::getTime() const
{
  return _time;
}

inline DiffusionMethod GrSimulation::getDiffusionMethod() const
{
  return _pDiffusion->getMethod();
}

inline void GrSimulation::setODESolverMethod(ODESolvers::ODEMethod odemethod)
{
  _odeSolver = odemethod;
}

inline int GrSimulation::getODESolverMethod() const
{
  return _odeSolver;
}

inline bool GrSimulation::getTnfrDynamics() const
{
  return _tnfrDynamics;
}

inline void GrSimulation::setTnfrDynamics(bool tnfrDynamics)
{
  _tnfrDynamics = tnfrDynamics;
}

inline bool GrSimulation::getNfkbDynamics() const
{
  return _nfkbDynamics;
}

inline void GrSimulation::setNfkbDynamics(bool nfkbDynamics)
{
  _nfkbDynamics = nfkbDynamics;
}

inline int GrSimulation::getTnfDepletionTimeStep() const
{
  return _tnfDepletionTimeStep;
}

inline void GrSimulation::setTnfDepletionTimeStep(int tnfDepletionTimeStep)
{
  _tnfDepletionTimeStep = tnfDepletionTimeStep;
}

inline bool GrSimulation::getIl10rDynamics() const
{
  return _il10rDynamics;
}

inline void GrSimulation::setIl10rDynamics(bool il10rDynamics)
{
  _il10rDynamics = il10rDynamics;
}

inline int GrSimulation::getIl10DepletionTimeStep() const
{
  return _il10DepletionTimeStep;
}

inline void GrSimulation::setIl10DepletionTimeStep(int il10DepletionTimeStep)
{
  _il10DepletionTimeStep = il10DepletionTimeStep;
}

inline bool GrSimulation::getTgammaTransition() const
{
  return _tgammatransition;
}

inline void GrSimulation::setTgammaTransition(bool tgammatransition)
{
  _tgammatransition = tgammatransition;
}
inline bool GrSimulation::getAdaptive() const
{
  return _adaptive;
}
inline void GrSimulation::setAdaptive(bool adaptive)
{
  _adaptive = adaptive;
}

inline const MacList& GrSimulation::getMacList() const
{
  return _macList;
}

inline const TgamList& GrSimulation::getTgamList() const
{
  return _tgamList;
}

inline const TcytList& GrSimulation::getTcytList() const
{
  return _tcytList;
}

inline const TregList& GrSimulation::getTregList() const
{
  return _tregList;
}

inline const GrGrid& GrSimulation::getGrid() const
{
  return _grid.getGrid();
}

inline GrGrid& GrSimulation::getGrid()
{
  return _grid.getGrid();
}

inline const Stats& GrSimulation::getStats() const
{
  return _stats;
}

inline Stats& GrSimulation::getStats()
{
  return _stats;
}

inline const Stats& GrSimulation::getStatsPrevious() const
{
  return _statsPrevious;
}

inline Stats& GrSimulation::getStatsPrevious()
{
  return _statsPrevious;
}

inline Mac* GrSimulation::createMac(int row, int col, int birthtime, Mac::State state, bool NFkB, bool stat1)
{
  _macList.push_back(new Mac(birthtime, row, col, state, 0, NFkB, stat1));
  Mac* pMac = _macList.back();

  assert_res(_grid.getGrid().addAgent(pMac, row, col));
  _stats.updateAgentStatistics(pMac);

  return pMac;
}

inline Tgam* GrSimulation::createTgam(int row, int col, int birthtime, Tgam::State state)
{
  _tgamList.push_back(new Tgam(birthtime, row, col, state));
  Tgam* pTgam = _tgamList.back();

  assert_res(_grid.getGrid().addAgent(pTgam, row, col));
  _stats.updateAgentStatistics(pTgam);

  return pTgam;
}

inline Tcyt* GrSimulation::createTcyt(int row, int col, int birthtime, Tcyt::State state)
{
  _tcytList.push_back(new Tcyt(birthtime, row, col, state));
  Tcyt* pTcyt = _tcytList.back();

  assert_res(_grid.getGrid().addAgent(pTcyt, row, col));
  _stats.updateAgentStatistics(pTcyt);

  return pTcyt;
}

inline Treg* GrSimulation::createTreg(int row, int col, int birthtime, Treg::State state)
{
  _tregList.push_back(new Treg(birthtime, row, col, state));
  Treg* pTreg = _tregList.back();

  assert_res(_grid.getGrid().addAgent(pTreg, row, col));
  _stats.updateAgentStatistics(pTreg);

  return pTreg;
}

inline void GrSimulation::checkTCellRecruitmentStart()
{
  if (_tcellRecruitmentBegun)
    {
      return;
    }

  if (_PARAM(PARAM_TCELL_MTB_RECRUITMENT_THRESHOLD) > 0)
    {
      Scalar totMtb = _stats.getTotExtMtb() + _stats.getTotIntMtb();
      if (totMtb > _PARAM(PARAM_TCELL_MTB_RECRUITMENT_THRESHOLD))
        {
          _tcellRecruitmentBegun = true;
        }
    }
  else
    {
      if (_time >= _PARAM(PARAM_TCELL_TIME_RECRUITMENT_ENABLED))
        {
          _tcellRecruitmentBegun = true;
        }
    }
}

inline void GrSimulation::timestepSync()
{
//    // Checks to see if the timesteps specified will work with a 10 minutes agent timestep
//    if (((SECONDS_PER_DAY)/TIME_STEPS_PER_DAY) % _PARAM(PARAM_GR_DT_DIFFUSION) != 0 || ((SECONDS_PER_DAY)/TIME_STEPS_PER_DAY) % _PARAM(PARAM_GR_DT_MOLECULAR) != 0 || _PARAM(PARAM_GR_DT_DIFFUSION) % _PARAM(PARAM_GR_DT_MOLECULAR) != 0)
//    {
//        std::cerr << "\nUnsupported Time Step:\n--Diffusion/Molecular time step must be able to sync with the 10 min agent time step\n" << std::endl;
//        throw std::runtime_error("*** ERROR: --Redefine the timesteps in the parameter file");
//    }
//    // Check the diffusion time step and see if it is larger than the heuristic for level of accuracy
//    if (_PARAM(PARAM_GR_DT_DIFFUSION) >= (2.0/(_PARAM(PARAM_GR_D_TNF) * ((1/pow(20e-4,2)) + (1/pow(20e-4,2))))))
//    {
//        std::cerr << "\nWarning:\n--The diffusion time step is higher than its approximate accuracy limit\n--Any time step higher than ~60 seconds will be a stable but more inaccurate solution\n" << std::endl;
//        throw std::runtime_error("*** ERROR: --Redefine the timesteps in the parameter file");
//
//    }
//    // Verify that the molecular time step is small than the diffusion time step and smaller than the approximate stability limit
//    if (_PARAM(PARAM_GR_DT_MOLECULAR) > MOLECULAR_ACCURACY || _PARAM(PARAM_GR_DT_MOLECULAR) > _PARAM(PARAM_GR_DT_DIFFUSION))
//    {
//        std::cerr << "\nUnsupported Molecular Time Step:\n--The solution is unstable at time steps greater than 30s\n--The molecular time step must be smaller than the diffusion time step\n" << std::endl;
//        throw std::runtime_error("*** ERROR: --Redefine the timesteps in the parameter file");
//    }

  double sync1int, sync2int, sync3int;
  double sync1frac = modf(((SECONDS_PER_DAY)/TIME_STEPS_PER_DAY)/_PARAM(PARAM_GR_DT_DIFFUSION), &sync1int);
  double sync2frac = modf(((SECONDS_PER_DAY)/TIME_STEPS_PER_DAY)/_PARAM(PARAM_GR_DT_MOLECULAR), &sync2int);
  double sync3frac = modf(_PARAM(PARAM_GR_DT_DIFFUSION)/_PARAM(PARAM_GR_DT_MOLECULAR), &sync3int);

  // Checks to see if the timesteps specified will work with a 10 minutes agent timestep
  if (sync1frac != 0.0 || sync2frac != 0.0 || sync3frac != 0.0)
    {
      std::cerr << "\nUnsupported Time Step:\n--Diffusion/Molecular time step must be able to sync with the 10 min agent time step\n" << std::endl;
      throw std::runtime_error("*** ERROR: --Redefine the timesteps in the parameter file");
    }
  // Check the diffusion time step and see if it is larger than the heuristic for level of accuracy
  if (_PARAM(PARAM_GR_DT_DIFFUSION) >= (2.0/(_PARAM(PARAM_GR_D_TNF) * ((1/pow(20e-4,2)) + (1/pow(20e-4,2))))))
    {
      std::cerr << "\nWarning:\n--The diffusion time step is higher than its approximate accuracy limit\n--Any time step higher than ~60 seconds will be a stable but more inaccurate solution\n" << std::endl;
      throw std::runtime_error("*** ERROR: --Redefine the timesteps in the parameter file");

    }
  // Verify that the molecular time step is small than the diffusion time step and smaller than the approximate stability limit
  if (_PARAM(PARAM_GR_DT_MOLECULAR) > MOLECULAR_ACCURACY || _PARAM(PARAM_GR_DT_MOLECULAR) > _PARAM(PARAM_GR_DT_DIFFUSION))
    {
      std::cerr << "\nUnsupported Molecular Time Step:\n--The solution is unstable at time steps greater than 30s\n--The molecular time step must be smaller than the diffusion time step\n" << std::endl;
      throw std::runtime_error("*** ERROR: --Redefine the timesteps in the parameter file");
    }
}


template<class Archive>
void GrSimulation::serialize(Archive& ar, const unsigned int version) {
  #ifdef SVN_VERSION
    if(version != SVN_VERSION % 256)
      { 
        std::cerr<< "Warning: Serialization version mismatch!  Trying to serialize "
                 << (SVN_VERSION/256)*256 + version << " but program is " << SVN_VERSION <<std::endl;
      }
  #endif
  /* Register the classes to be serialized (only needed for serializing through
   * abstract base pointer type, like AgentPtr */
  boost::mpl::for_each<AgentTypes, wrap<boost::mpl::placeholders::_1> >
    (serialization_class_register_t<Archive>(ar));

  ar & BOOST_SERIALIZATION_NVP(_time);
  ar & BOOST_SERIALIZATION_NVP(_areaThreshold);
  ar & BOOST_SERIALIZATION_NVP(_areaThresholdCellDensity);
  //ar & BOOST_SERIALIZATION_NVP(_diffusion);
  //ar & BOOST_SERIALIZATION_NVP(_recruitment);
	ar & BOOST_SERIALIZATION_NVP(_tnfrDynamics);
	ar & BOOST_SERIALIZATION_NVP(_nfkbDynamics);
  ar & BOOST_SERIALIZATION_NVP(_il10rDynamics);
  ar & BOOST_SERIALIZATION_NVP(_tgammatransition);
  ar & BOOST_SERIALIZATION_NVP(_odeSolver);
	ar & BOOST_SERIALIZATION_NVP(_tnfDepletionTimeStep);
  ar & BOOST_SERIALIZATION_NVP(_il10DepletionTimeStep);
  ar & BOOST_SERIALIZATION_NVP(_tcellRecruitmentBegun);
  ar & BOOST_SERIALIZATION_NVP(_numMolecularPerDiffusion);
  ar & BOOST_SERIALIZATION_NVP(_numDiffusionPerAgent);

  ar & BOOST_SERIALIZATION_NVP(_macList);
  ar & BOOST_SERIALIZATION_NVP(_tgamList);
  ar & BOOST_SERIALIZATION_NVP(_tregList);
  ar & BOOST_SERIALIZATION_NVP(_tcytList);

  ar & BOOST_SERIALIZATION_NVP(_grid);

  Agent::classSerialize(ar, version);

  ar & BOOST_SERIALIZATION_NVP(_statsPrevious);
  ar & BOOST_SERIALIZATION_NVP(_stats);
  ar & BOOST_SERIALIZATION_NVP(g_Rand);
}

#ifdef SVN_VERSION  //Add the svn version to the class
  BOOST_CLASS_VERSION(GrSimulation, SVN_VERSION % 256);
#endif

#endif /* GRSIMULATION_H */
