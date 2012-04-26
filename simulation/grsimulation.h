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

using namespace std;

class GrSimulation
{
private:
	static const std::string _ClassName;

	/*
	 * !!! If the data members change then the serialize and deserialize functions need to be updated !!!
	 */
	int _time;
	GrSimulationGrid _grid;
	MacList _macList;
	TgamList _tgamList;
	TcytList _tcytList;
	TregList _tregList;

	Stats _statsPrevious; // The stats at the start of a time step - from the end of the previous time step.
	Stats _stats;

	double _areaThreshold;
	double _areaThresholdCellDensity;
	GrDiffusion* _pDiffusion;
	TTest* _pTTest[NOUTCOMES];
	RecruitmentBase* _pRecruitment;
	bool _tnfrDynamics;
	bool _nfkbDynamics;
    bool _il10rDynamics;
    bool _tgammatransition;

    // ODE Solver
    int _odeSolver;
    
	// Inhibits tnf secretion if true and if not using tnfr dynamics.
	int _tnfDepletionTimeStep;
    int _il10DepletionTimeStep;

    // Whether or not TCell recruitment has been enabled.
    // Once enabled it stays enabled even if the criteria by which it became enabled changes.
    bool _tcellRecruitmentBegun;
    
    int _numMolecularPerDiffusion;
    int _numDiffusionPerAgent;

	void moveTcells();
	void moveMacrophages();
	void updateStates();
	void updateT_Test();
	void computeNextStates();
	void secreteFromMacrophages(bool tnfDepletion, bool il10Depletion, int mdt);
	void secreteFromTcells(bool tnfDepletion, bool il10Depletion, int mdt);
	void secreteFromCaseations(int mdt);
    void updateTNFDynamics(double dt);
    void updateIL10Dynamics(double dt);
    void updateTNFandIL10Dynamics(double dt);
    void updateNFkBandTNFandIL10Dynamics(double dt);
    void updateNFkBandTNFDynamics(double dt);
    void updateMolecularScaleRK4(double dt);
    void updateMolecularScaleRK2(double dt);
    void updateMolecularScaleFE(double dt);
    void updateMolecularScaleEPC(double dt);
	void adjustTNFDegradation(double dt);
    void adjustFauxDegradation(double dt);
	void growExtMtb();
	void shuffleCells();
	void checkTCellRecruitmentStart();

public:
	GrSimulation(const Pos& dim);
	~GrSimulation();
	void init(Scalar molecularTrackingRadius);
	void initMolecularTracking(Scalar molecularTrackingRadius);
	void solve();
	void performT_Test();
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
    void setODESolverMethod(int odemethod);
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
	double getAreaThreshold() const;
	void setAreaThreshold(double areaThreshold);
	double getAreaThresholdCellDensity() const;
	void setAreaThresholdCellDensity(double areaThreshold);
	OutcomeMethod getOutcomeMethod(int index) const;
	void getOutcomeParameters(int index, int& samplePeriod, int& testPeriod, double& alpha) const;
	void setOutcomeMethod(int index, OutcomeMethod method, double alpha, int testPeriod, int samplePeriod);
	bool getTCellRecruitmentBegun();
	void serialize(std::ostream& out) const;
	void deserialize(std::istream& in);
	void setRecruitment	(RecruitmentBase* pRecruitment);
	Mac* createMac(int row, int col, int birthtime, Mac::State state, bool NFkB, bool stat1);
	Tgam* createTgam(int row, int col, int birthtime, Tgam::State state);
	Tcyt* createTcyt(int row, int col, int birthtime, Tcyt::State state);
	Treg* createTreg(int row, int col, int birthtime, Treg::State state);

	static void convertSimTime(const int time, int& rDays, int& rHours, int& rMinutes);
    void setODEsize();
    void macODEsize();
    void tcellODEsize();
    void timestepSync();
};

inline bool GrSimulation::getTCellRecruitmentBegun()
{
	return _tcellRecruitmentBegun;
}

inline void GrSimulation::setRecruitment(RecruitmentBase* pRecruitment)
{
	_pRecruitment = pRecruitment;
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

inline void GrSimulation::setODESolverMethod(int odemethod)
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
	_macList.push_back(Mac(birthtime, row, col, state, 0, NFkB, stat1));
	Mac* pMac = &_macList.back();
	
	assert_res(_grid.getGrid().addAgent(pMac, row, col));
	_stats.updateAgentStatistics(pMac);

	return pMac;
}

inline Tgam* GrSimulation::createTgam(int row, int col, int birthtime, Tgam::State state)
{
	_tgamList.push_back(Tgam(birthtime, row, col, state));
	Tgam* pTgam = &_tgamList.back();
	
	assert_res(_grid.getGrid().addAgent(pTgam, row, col));
	_stats.updateAgentStatistics(pTgam);

	return pTgam;
}

inline Tcyt* GrSimulation::createTcyt(int row, int col, int birthtime, Tcyt::State state)
{
	_tcytList.push_back(Tcyt(birthtime, row, col, state));
	Tcyt* pTcyt = &_tcytList.back();
	
	assert_res(_grid.getGrid().addAgent(pTcyt, row, col));
	_stats.updateAgentStatistics(pTcyt);

	return pTcyt;
}

inline Treg* GrSimulation::createTreg(int row, int col, int birthtime, Treg::State state)
{
	_tregList.push_back(Treg(birthtime, row, col, state));
	Treg* pTreg = &_tregList.back();
	
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

inline void GrSimulation::setODEsize()
{
    macODEsize();
    tcellODEsize();
}

inline void GrSimulation::macODEsize()
{
    int vectorlength = 0;
    
    if (_nfkbDynamics)
    {
        if (_il10rDynamics) 
        {
            vectorlength = 41; // Number of differential equations for TNF, IL10, and NFKB
        }
        
        else
        {
            vectorlength = 38; // Number of differential equations for TNF and NFKB
        }
    }
    
    else if (_tnfrDynamics && _il10rDynamics)
    {
        vectorlength = 13; // Number of differential equations for TNF and IL10
    }
    
    else if (_tnfrDynamics && !_il10rDynamics)
    {
        vectorlength = 10; // Number of differential equations for TNF
    }
    
    else if (_il10rDynamics && !_tnfrDynamics)
    {
        vectorlength = 3; // Number of differential equations for IL10
    }
    
    Mac::setMacOdeSize(vectorlength);
    
}


inline void GrSimulation::tcellODEsize()
{
    int vectorlength = 0;
    
    if (_tnfrDynamics && _il10rDynamics)
    {
        vectorlength = 13; // Number of differential equations for TNF and IL10
    }
    
    else if (_tnfrDynamics && !_il10rDynamics)
    {
        vectorlength = 10; // Number of differential equations for TNF
    }
    
    else if (_il10rDynamics && !_tnfrDynamics)
    {
        vectorlength = 3; // Number of differential equations for IL10
    }
    
    Tcell::setTcellOdeSize(vectorlength);
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



#endif /* GRSIMULATION_H */
