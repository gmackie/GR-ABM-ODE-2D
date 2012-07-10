/*
 * grsimulation.cpp
 *
 *  Created on: 05-nov-2009
 *      Author: M. El-Kebir
 */

#include "grsimulation.h"
#include "float.h"
#include "grdiffusion.h"
#include "grdiffusionbtcs.h"
#include "grdiffusionwrongbtcs.h"
#include "grdiffusionadeswap.h"
#include "grdiffusionftcsswap.h"
#include "grdiffusionadeswap.h"
#include "areatest.h"
#include "mtbtest.h"
#include "recruitmentprob.h"
#include "recruitmentlnodepure.h"
#include "serialization.h"
#include <boost/date_time/posix_time/posix_time.hpp>

using namespace std;

const std::string GrSimulation::_ClassName = "GrSimulation";

GrSimulation::GrSimulation(const Pos& dim)
	: _time(0)
	, _grid(dim)
	, _macList()
	, _tgamList()
	, _tcytList()
	, _tregList()
	, _statsPrevious()
	, _stats()
	, _areaThreshold(0.5f)
	, _areaThresholdCellDensity(0.5f)
	, _pDiffusion(new GrDiffusionADE_Swap())
	, _pTTest()
	, _pRecruitment(NULL)
	, _tnfrDynamics(false)
	, _nfkbDynamics(false)
  , _il10rDynamics(false)
  , _tgammatransition(false)
  , _odeSolver(ODESolvers::NMethods)
  , _tnfDepletionTimeStep(-1)
  , _il10DepletionTimeStep(-1)
	, _tcellRecruitmentBegun(false)
  , _numMolecularPerDiffusion(0)
  , _numDiffusionPerAgent(0)
{
	for (int i = 0; i < NOUTCOMES; i++)
    {
		_pTTest[i] = NULL;
    }

    _numMolecularPerDiffusion = _PARAM(PARAM_GR_DT_DIFFUSION)/_PARAM(PARAM_GR_DT_MOLECULAR); // Number of molecular iterations per diffusion iteration
    _numDiffusionPerAgent = (AGENT_TIME_STEP / _PARAM(PARAM_GR_DT_DIFFUSION)); // Number of diffusion iterations per agent iteration
}

GrSimulation* GrSimulation::clone() const {
  GrSimulation* ret = new GrSimulation(*this);
  for(MacList::iterator begin = ret->_macList.begin(); begin!=ret->_macList.end(); begin++)
    (*begin) = static_cast<Mac*>((*begin)->clone());
  for(TgamList::iterator begin = ret->_tgamList.begin(); begin!=ret->_tgamList.end(); begin++)
    (*begin) = static_cast<Tgam*>((*begin)->clone());
  for(TregList::iterator begin = ret->_tregList.begin(); begin!=ret->_tregList.end(); begin++)
    (*begin) = static_cast<Treg*>((*begin)->clone());
  for(TcytList::iterator begin = ret->_tcytList.begin(); begin!=ret->_tcytList.end(); begin++)
    (*begin) = static_cast<Tcyt*>((*begin)->clone());
  if(ret->_pDiffusion)
    ret->_pDiffusion = ret->_pDiffusion->clone();
  if(ret->_pRecruitment)
    ret->_pRecruitment = ret->_pRecruitment->clone();
  for(size_t i=0;i<NOUTCOMES;i++)
    if(ret->_pTTest[i])
      ret->_pTTest[i] = ret->_pTTest[i]->clone();
  return ret;
}


template<typename T>
inline void clearPtrList(T& container) {
  typename T::iterator end = container.end();
  for(typename T::iterator begin = container.begin(); begin!=end; begin++)
    delete *begin;
  container.clear();
}

GrSimulation::~GrSimulation()
{
	delete _pDiffusion;
	for (int i = 0; i < NOUTCOMES; i++)
		delete _pTTest[i];
  clearPtrList(_macList);
  clearPtrList(_tgamList);
  clearPtrList(_tregList);
  clearPtrList(_tcytList);
}

void GrSimulation::serialize(std::ostream& out) const
{
	assert(out.good());

	/*** Set initial fmt flags ***/
	out.setf(std::ios::scientific, std::ios::floatfield);
	out.precision(Scalar_PRECISION);

	Serialization::writeHeader(out, GrSimulation::_ClassName);

	// Model version
	out << GR_VERSION << std::endl;

	// serialize time
	out << _time << std::endl;

	// serialize area threshold
	out << _areaThreshold << std::endl;
	out << _areaThresholdCellDensity << std::endl;

	// serialize diffusion method
	int intVal = (int) _pDiffusion->getMethod();
	out << intVal << std::endl;

	// serialize recruitment method and object.
	int intRecrMtehod = (int) _pRecruitment->getMethod();
	out << intRecrMtehod << std::endl;
	_pRecruitment->serialize(out);

	out << _tnfrDynamics << std::endl;
	out << _nfkbDynamics << std::endl;
    out << _il10rDynamics << std::endl;
    out << _tgammatransition << std::endl;
    out << (int)(_odeSolver) << std::endl;
	out << _tnfDepletionTimeStep << std::endl;
    out << _il10DepletionTimeStep << std::endl;
    out << _tcellRecruitmentBegun <<std::endl;
    out << _numMolecularPerDiffusion << std::endl;
    out << _numDiffusionPerAgent << std::endl;

	// serialize grid
	_grid.serialize(out);

	// Serialize the agent class.
	Agent::classSerialize(out);

	// serialize macs
	out << _macList.size() << std::endl;
	for (MacList::const_iterator it = _macList.begin(); it != _macList.end(); it++)
	{
		int agentIndex = getGrid().agentIndex(*it);

		if (agentIndex < 0 || agentIndex >= (int) GrGrid::MAX_AGENTS_PER_CELL)
		{
			cerr << "Mac serialization error, invalid agentIndex: " << agentIndex
				 << " Mac address: " << (*it)
				 << " Mac position: " << (*it)->getPosition()
				 << endl;
			exit(1);
		}

		out << getGrid().agentIndex((*it)) << endl;
		(*it)->serialize(out);
	}

	// serialize tgam cells
	out << _tgamList.size() << std::endl;
	for (TgamList::const_iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
	{
		int agentIndex = getGrid().agentIndex(*it);

		if (agentIndex < 0 || agentIndex >= (int) GrGrid::MAX_AGENTS_PER_CELL)
		{
			cerr << "Tgam serialization error, invalid agentIndex: " << agentIndex
				 << " Tgam address: " << (*it)
				 << " Tgam position: " << (*it)->getPosition()
				 << endl;
			exit(1);
		}

		out << getGrid().agentIndex(*it) << endl;
		(*it)->serialize(out);
	}

	// serialize tcyt cells
	out << _tcytList.size() << std::endl;
	for (TcytList::const_iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
		int agentIndex = getGrid().agentIndex((*it));

		if (agentIndex < 0 || agentIndex >= (int) GrGrid::MAX_AGENTS_PER_CELL)
		{
			cerr << "Tcyt serialization error, invalid agentIndex: " << agentIndex
				 << " Tcyt address: " << &(*it)
				 << " Tcyt position: " << (*it)->getPosition()
				 << endl;
			exit(1);
		}

		out << getGrid().agentIndex((*it)) << endl;
		(*it)->serialize(out);
	}

	// serialize treg cells
	out << _tregList.size() << std::endl;
	for (TregList::const_iterator it = _tregList.begin(); it != _tregList.end(); it++)
	{
		int agentIndex = getGrid().agentIndex((*it));

		if (agentIndex < 0 || agentIndex >= (int) GrGrid::MAX_AGENTS_PER_CELL)
		{
			cerr << "Treg serialization error, invalid agentIndex: " << agentIndex
				 << " Treg address: " << &(*it)
				 << " Treg position: " << (*it)->getPosition()
				 << endl;
			exit(1);
		}

		out << getGrid().agentIndex((*it)) << endl;
		(*it)->serialize(out);
	}

	// serialize statistics
	_statsPrevious.serialize(out);
	_stats.serialize(out);

	// serialize random number generator
	// Do this last so that when de-serializing any random number generation performed as part
	// de-serialization won't affect the simulation state.
	g_Rand.serialize(out);

	Serialization::writeFooter(out, GrSimulation::_ClassName);
}

void GrSimulation::deserialize(std::istream& in)
{

  timestepSync();

  setODEsize();

	assert(in.good());

	if (!Serialization::readHeader(in, GrSimulation::_ClassName))
	{
		exit(1);
	}

	std::string svnVersion;
	in >> svnVersion;

	// deserialize time
	in >> _time;

	// deserialize area threshold
	in >> _areaThreshold;
	in >> _areaThresholdCellDensity;

	// deserialize diffusion method
	int intVal;
	in >> intVal;
	setDiffusionMethod((DiffusionMethod) intVal);

	// deserialize recruitment method and object.
	int intRecrMtehod;
	in >> intRecrMtehod;
	RecruitmentMethod recrmethod = (RecruitmentMethod) intRecrMtehod;
	deserializeRecruitmentMethod(recrmethod, in);
  assert(in.good());
	_pRecruitment->deserialize(in); //LnODE recruitment will deserialize twice.
  assert(in.good());

  assert(in.good());
	in >> _tnfrDynamics;
	in >> _nfkbDynamics;
    in >> _il10rDynamics;
    in >> _tgammatransition;
    in >> (int&)(_odeSolver);
	in >> _tnfDepletionTimeStep;
    in >> _il10DepletionTimeStep;
    in >> _tcellRecruitmentBegun;
    in >> _numMolecularPerDiffusion;
    in >> _numDiffusionPerAgent;
  assert(in.good());

	// deserialize grid
	_grid.deserialize(in);

	// Deserialize the agent class.
	Agent::classDeserialize(in);

	// Since the cell lists are cleared here, the grid must not contain
	// any cell pointers at this point. If it does they are dangling pointers
	// that will be referenced during cell deserialization, when determining
	// if the grid compartment for the deserialized cell has room for the
	// deserialized cell. This will cause undefined (bad) behavior - seg faults, etc.

	// deserialize macs
	_macList.clear();
	in >> intVal;
	for (int i = 0; i < intVal; i++)
	{
        // Get the index for the grid compartment agent array.
        int index;
        in >> index;

		// create a dummy mac
		_macList.push_back(new Mac());
		Mac* pMac = _macList.back();

		// update attributes of dummy mac and add to grid
		pMac->deserialize(in);

		if (index < 0 || index >= (int) GrGrid::MAX_AGENTS_PER_CELL)
		{
			cerr << "Error deserializing mac " << i << " at grid position " << pMac->getPosition() << ": invalid grid agent array index " << index << endl;
			exit(1);
		}

		bool addedAgent = _grid.getGrid().addAgent(pMac, index);
		if (!addedAgent)
		{
			cerr << "Error deserializing mac " << i << " at grid position " << pMac->getPosition() << ": addAgent failed."<< endl;
			exit(1);
		}
	}

	// deserialize tgam cells
	_tgamList.clear();
	in >> intVal;
	for (int i = 0; i < intVal; i++)
	{
        // Get the index for the grid compartment agent array.
        int index;
        in >> index;

		// create a dummy tgam cell
		_tgamList.push_back(new Tgam());
		Tgam* pTgam = _tgamList.back();

		if (index < 0 || index >= (int) GrGrid::MAX_AGENTS_PER_CELL)
		{
			cerr << "Error deserializing tgam " << i << " at grid position " << pTgam->getPosition() << ": invalid grid agent array index " << index << endl;
			exit(1);
		}

		// update attributes of dummy tgam and add to grid
		pTgam->deserialize(in);
		bool addedAgent = _grid.getGrid().addAgent(pTgam, index);
		if (!addedAgent)
		{
			cerr << "Error deserializing tgam " << i << " at grid position " << pTgam->getPosition() << ": addAgent failed."<< endl;
			exit(1);
		}
	}

	// deserialize tcyt cells
	_tcytList.clear();
	in >> intVal;
	for (int i = 0; i < intVal; i++)
	{
        // Get the index for the grid compartment agent array.
        int index;
        in >> index;

		// create a dummy tcyt cell
		_tcytList.push_back(new Tcyt());
		Tcyt* pTcyt = _tcytList.back();

		if (index < 0 || index >= (int) GrGrid::MAX_AGENTS_PER_CELL)
		{
			cerr << "Error deserializing tcyt " << i << " at grid position " << pTcyt->getPosition() << ": invalid grid agent array index " << index << endl;
			exit(1);
		}

		// update attributes of dummy tcyt and add to grid
		pTcyt->deserialize(in);
		bool addedAgent = _grid.getGrid().addAgent(pTcyt, index);
		if (!addedAgent)
		{
			cerr << "Error deserializing tcyt " << i << " at grid position " << pTcyt->getPosition() << ": addAgent failed."<< endl;
			exit(1);
		}

	}

	// deserialize treg cells
	_tregList.clear();
	in >> intVal;
	for (int i = 0; i < intVal; i++)
	{
        // Get the index for the grid compartment agent array.
        int index;
        in >> index;

		// create a dummy treg cell
		_tregList.push_back(new Treg());
		Treg* pTreg = _tregList.back();

		if (index < 0 || index >= (int) GrGrid::MAX_AGENTS_PER_CELL)
		{
			cerr << "Error deserializing treg " << i << " at grid position " << pTreg->getPosition() << ": invalid grid agent array index " << index << endl;
			exit(1);
		}

		// update attributes of dummy mac and add to grid
		pTreg->deserialize(in);
		bool addedAgent = _grid.getGrid().addAgent(pTreg, index);
		if (!addedAgent)
		{
			cerr << "Error deserializing treg " << i << " at grid position " << pTreg->getPosition() << ":addAgent failed."<< endl;
			exit(1);
		}

	}

	// deserialize statistics
	_statsPrevious.serialize(in);
	_stats.serialize(in);

	// deserialize random number generator
	g_Rand.deserialize(in);

	if (!Serialization::readFooter(in, GrSimulation::_ClassName))
	{
		exit(1);
	}
}

void GrSimulation::init()
{
	// Before initializing anything check to see if the time step criteria are met
    timestepSync();
    
    // Used to initialize the valarray length for solving ODEs using RK4
    // Must be run before any agents are created so they have the correct
    // valarray length
    setODEsize();

    // initialize the sources
	_grid.initSources();

	PosVector initMacs = Params::getInstance()->getInitialMacs();
	PosVector initExtMtb = Params::getInstance()->getInitialExtMtb();
	const int maxMacAge = _PARAM(PARAM_MAC_AGE);

	// Place initial infected macrophages on the grid
	for (PosVector::const_iterator it = initMacs.begin(); it != initMacs.end(); it++)
	{
		Mac* pMac = createMac(it->x, it->y, g_Rand.getInt(-1, -maxMacAge), Mac::MAC_INFECTED, false, false);
		
		if (_nfkbDynamics)
		{
			// initialize NF-kB signaling from steady-state (12 hours equilibrium)
			for (int i = 0; i < 7200*_PARAM(PARAM_GR_NF_KB_TIME_COEFF); ++i)
				pMac->solveNFkBODEsEquilibrium(6.00/_PARAM(PARAM_GR_NF_KB_TIME_COEFF));
		}

		pMac->setIntMtb(1);
		++_statsPrevious.getTotIntMtb();
		++_stats.getTotIntMtb();
	}

	// Place initial extracellular bacteria on the grid
	for (PosVector::const_iterator it = initExtMtb.begin(); it != initExtMtb.end(); it++)
	{
		_grid.getGrid().extMTB(*it) += 1;
		++_statsPrevious.getTotExtMtb();
		++_stats.getTotExtMtb();
	}

	// Add a fixed number (PARAM_MAC_INIT_NUMBER) of resting
	// macrophages to the grid
	int count = _PARAM(PARAM_MAC_INIT_NUMBER);

	int nrGridCompartments = (_grid.getSize());
	if (count > nrGridCompartments)
	{
		std::cerr << "The number of initial resting macrophages to place on the grid, " << count << ", is > the number of grid compartments, " << nrGridCompartments << "." << std::endl;
		exit(1);
	}

	while (count > 0)
	{
		int row = g_Rand.getInt(_grid.getRange().x);
		int col = g_Rand.getInt(_grid.getRange().y);
		if (!_grid.getGrid().hasAgentType(MAC, Pos(row, col)))
		{
			Mac* newMac = createMac(row, col, g_Rand.getInt(-1, -maxMacAge), Mac::MAC_RESTING, false, false);
			if (_nfkbDynamics)
			{
				// initialize NF-kB signaling from steady-state (12 hours equilibrium)
				for (int i = 0; i < 7200*_PARAM(PARAM_GR_NF_KB_TIME_COEFF); ++i)
					newMac->solveNFkBODEsEquilibrium(6.00/_PARAM(PARAM_GR_NF_KB_TIME_COEFF));
			}
			count--;
		}
	}

}

// Mark each cell within the molecular tracking radius to be tracked.
// Clear a cell's tracking status otherwise, so that the cells to track for a loaded
// state doesn't include those that were marked to be tracked when the state was saved.
// This allows defining the cells to be tracked based on the simulation state at the
// time the saved state was loaded, rather than based in part on the simulation state
// at time 0, when the simulation was initialized, and which would be carried forward to
// the saved and then loaded state.
void GrSimulation::initMolecularTracking(Scalar molecularTrackingRadius)
{
	if (molecularTrackingRadius == 0.0)
	{
		return;
	}

	const Pos center = _grid.getCenter();

	for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
	{
		if ((*it)->getPosition().distance(center) <= molecularTrackingRadius)
		{
			(*it)->setTrackMolecularDynamics(true);
		}
		else
		{
			(*it)->setTrackMolecularDynamics(false);
		}

	}

	for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
	{
		if ((*it)->getPosition().distance(center) <= molecularTrackingRadius)
		{
			(*it)->setTrackMolecularDynamics(true);
		}
		else
		{
			(*it)->setTrackMolecularDynamics(false);
		}
	}

	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
		if ((*it)->getPosition().distance(center) <= molecularTrackingRadius)
		{
			(*it)->setTrackMolecularDynamics(true);
		}
		else
		{
			(*it)->setTrackMolecularDynamics(false);
		}
	}

	for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
	{
		if ((*it)->getPosition().distance(center) <= molecularTrackingRadius)
		{
			(*it)->setTrackMolecularDynamics(true);
		}
		else
		{
			(*it)->setTrackMolecularDynamics(false);
		}
	}
}
void GrSimulation::initMolecularTracking(const std::vector<size_t>& ids) {

	for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
    for(typeof(ids.begin()) i=ids.begin(); i != ids.end(); i++)
      if(*i == (*it)->getID())
        (*it)->setTrackMolecularDynamics(true);
	for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
    for(typeof(ids.begin()) i=ids.begin(); i != ids.end(); i++)
      if(*i == (*it)->getID())
        (*it)->setTrackMolecularDynamics(true);
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
    for(typeof(ids.begin()) i=ids.begin(); i != ids.end(); i++)
      if(*i == (*it)->getID())
        (*it)->setTrackMolecularDynamics(true);
	for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
    for(typeof(ids.begin()) i=ids.begin(); i != ids.end(); i++)
      if(*i == (*it)->getID())
        (*it)->setTrackMolecularDynamics(true);
}

void GrSimulation::solve()
{
  using namespace boost::posix_time;
	// here we perform a timestep, which is 10 minutes
	_time++;
  ptime startWallTime(microsec_clock::local_time());

	#if 0
	//DBG
	// This is useful for determining whether or not 2 versions of the code are
	// using the same random number sequence, if they are using the same seed
	// and run on the same system.
	g_Rand.test(_time);
	//DBG
	#endif

	// Ensure that both grids have the same state.
	_grid.updateNextGrid();

	// reset statistics
	_stats.reset();

	bool tnfDepletion = false;
    bool il10Depletion = false;

	if (_tnfDepletionTimeStep >= 0 && _time >= _tnfDepletionTimeStep)
	{
		tnfDepletion = true;
		_tnfrDynamics = false;
	}

    if (_il10DepletionTimeStep >= 0 && _time >= _il10DepletionTimeStep)
	{
		il10Depletion = true;
		_il10rDynamics = false;
	}
    
    // Calculate Diffusion and Molecular events for a 10 min timestep
    // Agent time step is 600 s (10 min)
	
	for (int DiffStep = 0; DiffStep < _numDiffusionPerAgent; DiffStep++) 
	{
		_pDiffusion->diffuse(_grid);
                
#if 0
        if (_odeSolver == 4) 
        {
            secreteFromMacrophages(tnfDepletion, il10Depletion, _PARAM(PARAM_GR_DT_DIFFUSION));
            secreteFromTcells(tnfDepletion, il10Depletion, _PARAM(PARAM_GR_DT_DIFFUSION));
            secreteFromCaseations(_PARAM(PARAM_GR_DT_DIFFUSION));
            
            if (_nfkbDynamics || _tnfrDynamics || _il10rDynamics)
            {
              assert(0 && "NOT IMPLEMENTED");
                //updateMolecularScaleAdaptive(_PARAM(PARAM_GR_DT_MOLECULAR));
            }
            
            else if (!_il10rDynamics || !_tnfrDynamics)
            {
                adjustFauxDegradation(_PARAM(PARAM_GR_DT_DIFFUSION));
                //adjustTNFDegradation(_PARAM(PARAM_GR_DT_MOLECULAR));
            }
        }
        
        else
        {
#endif
            for (int MolStep = 0; MolStep < _numMolecularPerDiffusion; MolStep++)
            {
                secreteFromMacrophages(tnfDepletion, il10Depletion, _PARAM(PARAM_GR_DT_MOLECULAR));
                secreteFromTcells(tnfDepletion, il10Depletion, _PARAM(PARAM_GR_DT_MOLECULAR));
                secreteFromCaseations(_PARAM(PARAM_GR_DT_MOLECULAR));

                if (_nfkbDynamics || _tnfrDynamics || _il10rDynamics)
                  if(_adaptive)
                    solveMolecularScaleAdaptive(_PARAM(PARAM_GR_DT_MOLECULAR));
                  else
                    solveMolecularScale(_PARAM(PARAM_GR_DT_MOLECULAR));
                
#if 0
                if (_nfkbDynamics || _tnfrDynamics || _il10rDynamics) {
                    switch (_odeSolver) {
                        case 0:
                            updateMolecularScaleFE(_PARAM(PARAM_GR_DT_MOLECULAR));
                            break;
                        case 1:
                            updateMolecularScaleRK4(_PARAM(PARAM_GR_DT_MOLECULAR));
                            break;
                        case 2:
                            updateMolecularScaleRK2(_PARAM(PARAM_GR_DT_MOLECULAR));
                            break;
                        case 3:
                            updateMolecularScaleEPC(_PARAM(PARAM_GR_DT_MOLECULAR));
                            break;
                        default:
                            std::cerr<<"Error: Solver Switch Value Unsupported"<<std::endl;
                            exit(1);
                            break;
                    }
                }
#endif
                //                if (_nfkbDynamics)
                //                {
                //                    if (_il10rDynamics) {
                //                        updateNFkBandTNFandIL10Dynamics(_PARAM(PARAM_GR_DT_MOLECULAR));
                //                    }
                //                    else
                //                    {
                //                        updateNFkBandTNFDynamics(_PARAM(PARAM_GR_DT_MOLECULAR));
                //                    }
                //                }
                //                else if (_tnfrDynamics && _il10rDynamics)
                //                {
                //                    updateTNFandIL10Dynamics(_PARAM(PARAM_GR_DT_MOLECULAR));
                //                }
                //                
                //                else if (_tnfrDynamics && !_il10rDynamics)
                //                {
                //                    updateTNFDynamics(_PARAM(PARAM_GR_DT_MOLECULAR));
                //                }
                //            
                //                else if (_il10rDynamics && !_tnfrDynamics)
                //                {
                //                    updateIL10Dynamics(_PARAM(PARAM_GR_DT_MOLECULAR));
                //                }
                
                else if (!_il10rDynamics || !_tnfrDynamics)
                {
                    adjustFauxDegradation(_PARAM(PARAM_GR_DT_MOLECULAR));
                    //adjustTNFDegradation(_PARAM(PARAM_GR_DT_MOLECULAR));
                }
        }  

    }
	
	// move macrophages
	moveMacrophages();

	// move T cells every 10 minutes
	moveTcells();
	
	// Shuffle the sources and recruit agents from vascular sources every 10 minutes
	getGrid().shuffleSources();
	_pRecruitment->recruit(*this);

	// reset statistics
	_stats.resetAgentStats();

	// compute next state every 10 minutes
	computeNextStates();

	// update extracellular Mtb
	growExtMtb();

	// update states and remove dead agents from lists and grid
	updateStates();

	// This must be after growExtMtb and updateStates, since updateStates updates the stats with intMtb count
	// and growExtMtb updates the stats with extMtb counts.
	// The mtb counts are set to zero at the start of function solve above.
	checkTCellRecruitmentStart();
	
	// shuffle order of cells every hour for randomization in movement
	if (_time % 6 == 0)
	{
		shuffleCells();
	}

	// perform t-test
	updateT_Test();

	// Copy statistics, so we can use them on the next time step.
	_statsPrevious = _stats;

  ptime endWallTime(microsec_clock::local_time());
  time_duration dur = (endWallTime - startWallTime);
  _stats.setTimePerStep(dur.total_microseconds()*1.0e-6);
}

void GrSimulation::updateStates()
{
	// update states and remove dead macrophages from lists
	for (MacList::iterator it = _macList.begin(); it != _macList.end();)
	{
		Mac& mac = **it;
		mac.updateState();
		_stats.updateAgentStatistics(&mac);

		if (mac.isDead())
		{
			assert_res(_grid.getGrid().removeAgent(&mac));
      delete *it;
			it = _macList.erase(it);
		}
		else
		{
			it++;
		}
	}

	// update states and remove dead Tgams from lists
	for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end();)
	{
		Tgam& tgam = **it;
		tgam.updateState();
		_stats.updateAgentStatistics(&tgam);

		if (tgam.isDead())
		{
			assert_res(_grid.getGrid().removeAgent(&tgam));
      delete *it;
			it = _tgamList.erase(it);
		}
		else
		{
			it++;
		}
	}

	// update states and remove dead Tcyts from lists
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end();)
	{
		Tcyt& tcyt = **it;
		tcyt.updateState();
		_stats.updateAgentStatistics(&tcyt);

		if (tcyt.isDead())
		{
			assert_res(_grid.getGrid().removeAgent(&tcyt));
      delete *it;
			it = _tcytList.erase(it);
		}
		else
		{
			it++;
		}
	}

	// update states and remove dead Tregs from lists
	for (TregList::iterator it = _tregList.begin(); it != _tregList.end();)
	{
		Treg& treg = **it;
		treg.updateState();
		_stats.updateAgentStatistics(&treg);

		if ((*it)->isDead())
		{
			assert_res(_grid.getGrid().removeAgent(&treg));
      delete *it;
			it = _tregList.erase(it);
		}
		else
		{
			it++;
		}
	}
}

void GrSimulation::computeNextStates()
{
	for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
	{
		(*it)->computeNextState(_time, _grid.getGrid(), _stats, _tnfrDynamics, _nfkbDynamics, _il10rDynamics, _tgammatransition);
	}
	for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
	{
		(*it)->computeNextState(_time, _grid.getGrid(), _stats, _tnfrDynamics, _nfkbDynamics, _il10rDynamics, _tgammatransition);
	}
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
		(*it)->computeNextState(_time, _grid.getGrid(), _stats, _tnfrDynamics, _nfkbDynamics, _il10rDynamics, _tgammatransition);
	}
	for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
	{
		(*it)->computeNextState(_time, _grid.getGrid(), _stats, _tnfrDynamics, _nfkbDynamics, _il10rDynamics, _tgammatransition);
	} 
}

void GrSimulation::secreteFromMacrophages(bool tnfDepletion, bool il10Depletion, int mdt)
{
	for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
	{
		(*it)->secrete(_grid.getGrid(), _tnfrDynamics, _nfkbDynamics, tnfDepletion, _il10rDynamics, il10Depletion, mdt);
	}
}

void GrSimulation::secreteFromTcells(bool tnfDepletion, bool il10Depletion, int mdt)
{
	for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
	{
		(*it)->secrete(_grid.getGrid(), _tnfrDynamics, _nfkbDynamics, tnfDepletion, _il10rDynamics, il10Depletion, mdt);
	}
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
		(*it)->secrete(_grid.getGrid(), _tnfrDynamics, _nfkbDynamics, tnfDepletion, _il10rDynamics, il10Depletion, mdt);
	}
    for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
	{
		(*it)->secrete(_grid.getGrid(), _tnfrDynamics, _nfkbDynamics, tnfDepletion, _il10rDynamics, il10Depletion, mdt);
	}
}

void GrSimulation::secreteFromCaseations(int mdt)
{
	// secrete chemokines from caseated compartments only if infection is not cleared
  Pos p;
  GrGrid& g = _grid.getGrid();
  const Pos& dim = g.getRange();
	for (p.x = 0; p.x < dim.x; p.x++)
	{
		for (p.y = 0; p.y < dim.y; p.y++)
		{
			if (g.isCaseated(p))
			{
				g.incCCL2(p, (0.25 * _PARAM(PARAM_MAC_SEC_RATE_CCL2) * mdt));
				g.incCCL5(p, (0.25 * _PARAM(PARAM_MAC_SEC_RATE_CCL5) * mdt));
				g.incCXCL9(p,  (0.25 * _PARAM(PARAM_MAC_SEC_RATE_CXCL9) * mdt));
				g.incmacAttractant(p, (_PARAM(PARAM_GR_SEC_RATE_ATTRACTANT) * mdt));
			}
		}
	}
}

void GrSimulation::solveMolecularScale(double dt)
{
  for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
        (*it)->solveMolecularScale(_grid.getGrid(), 0, dt, _odeSolver);   //t is not correct, should be getTime()*600+SEC_SINCE_LAST
	for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
        (*it)->solveMolecularScale(_grid.getGrid(), 0, dt, _odeSolver);
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
        (*it)->solveMolecularScale(_grid.getGrid(), 0, dt, _odeSolver);
	for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
        (*it)->solveMolecularScale(_grid.getGrid(), 0, dt, _odeSolver);
  
}
void GrSimulation::solveMolecularScaleAdaptive(double dt)
{
  for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
    (*it)->solveMolecularScaleAdaptive(_grid.getGrid(), 0, dt, _odeSolver);   //t is not correct, should be getTime()*600+SEC_SINCE_LAST
	for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
    if(!(*it)->getODEstatus())
      (*it)->solveMolecularScaleAdaptive(_grid.getGrid(), 0, dt, _odeSolver);
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
    if(!(*it)->getODEstatus())
      (*it)->solveMolecularScaleAdaptive(_grid.getGrid(), 0, dt, _odeSolver);
	for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
    if(!(*it)->getODEstatus())
      (*it)->solveMolecularScaleAdaptive(_grid.getGrid(), 0, dt, _odeSolver);

  for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
    (*it)->setODEstatus(false);
	for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
    (*it)->setODEstatus(false);
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
    (*it)->setODEstatus(false);
	for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
    (*it)->setODEstatus(false);
}

#if 0
void GrSimulation::updateMolecularScaleAdaptive(double dt)
{
    for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
	{
        (*it)->solveMolecularScaleRKadaptive(_grid.getGrid(), dt, _nfkbDynamics, _tnfrDynamics, _il10rDynamics, _PARAM(PARAM_GR_DT_DIFFUSION));
	}
	for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
	{
        (*it)->solveMolecularScaleRKadaptive(_grid.getGrid(), dt, _nfkbDynamics, _tnfrDynamics, _il10rDynamics, _PARAM(PARAM_GR_DT_DIFFUSION));
	}
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
        (*it)->solveMolecularScaleRKadaptive(_grid.getGrid(), dt, _nfkbDynamics, _tnfrDynamics, _il10rDynamics, _PARAM(PARAM_GR_DT_DIFFUSION));
	}
	for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
	{
        (*it)->solveMolecularScaleRKadaptive(_grid.getGrid(), dt, _nfkbDynamics, _tnfrDynamics, _il10rDynamics, _PARAM(PARAM_GR_DT_DIFFUSION));
	}
}


void GrSimulation::updateMolecularScaleRK4(double dt)
{
  for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
	{
        (*it)->solveMolecularScaleRK4(_grid.getGrid(), dt, _nfkbDynamics, _tnfrDynamics, _il10rDynamics);
	}
	for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
	{
        (*it)->solveMolecularScaleRK4(_grid.getGrid(), dt, _nfkbDynamics, _tnfrDynamics, _il10rDynamics);
	}
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
        (*it)->solveMolecularScaleRK4(_grid.getGrid(), dt, _nfkbDynamics, _tnfrDynamics, _il10rDynamics);
	}
	for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
	{
        (*it)->solveMolecularScaleRK4(_grid.getGrid(), dt, _nfkbDynamics, _tnfrDynamics, _il10rDynamics);
	}
}

void GrSimulation::updateMolecularScaleRK2(double dt)
{
    for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
	{
        (*it)->solveMolecularScaleRK2(_grid.getGrid(), dt, _nfkbDynamics, _tnfrDynamics, _il10rDynamics);
	}
	for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
	{
        (*it)->solveMolecularScaleRK2(_grid.getGrid(), dt, _nfkbDynamics, _tnfrDynamics, _il10rDynamics);
	}
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
        (*it)->solveMolecularScaleRK2(_grid.getGrid(), dt, _nfkbDynamics, _tnfrDynamics, _il10rDynamics);
	}
	for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
	{
        (*it)->solveMolecularScaleRK2(_grid.getGrid(), dt, _nfkbDynamics, _tnfrDynamics, _il10rDynamics);
	}
}

void GrSimulation::updateMolecularScaleFE(double dt)
{
    for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
	{
        (*it)->solveMolecularScaleFE(_grid.getGrid(), dt, _nfkbDynamics, _tnfrDynamics, _il10rDynamics);
	}
	for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
	{
        (*it)->solveMolecularScaleFE(_grid.getGrid(), dt, _nfkbDynamics, _tnfrDynamics, _il10rDynamics);
	}
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
        (*it)->solveMolecularScaleFE(_grid.getGrid(), dt, _nfkbDynamics, _tnfrDynamics, _il10rDynamics);
	}
	for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
	{
        (*it)->solveMolecularScaleFE(_grid.getGrid(), dt, _nfkbDynamics, _tnfrDynamics, _il10rDynamics);
	}
}

void GrSimulation::updateMolecularScaleEPC(double dt)
{
    for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
	{
        (*it)->solveMolecularScaleFE(_grid.getGrid(), dt, _nfkbDynamics, _tnfrDynamics, _il10rDynamics);
	}
	for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
	{
        (*it)->solveMolecularScaleFE(_grid.getGrid(), dt, _nfkbDynamics, _tnfrDynamics, _il10rDynamics);
	}
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
        (*it)->solveMolecularScaleFE(_grid.getGrid(), dt, _nfkbDynamics, _tnfrDynamics, _il10rDynamics);
	}
	for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
	{
        (*it)->solveMolecularScaleFE(_grid.getGrid(), dt, _nfkbDynamics, _tnfrDynamics, _il10rDynamics);
	}
}

void GrSimulation::updateTNFDynamics(double dt)
{
	for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
	{
		(*it)->solveTNF(_grid.getGrid(), dt);
	}
	for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
	{
		(*it)->solveTNF(_grid.getGrid(), dt);
	}
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
		(*it)->solveTNF(_grid.getGrid(), dt);
	}
	for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
	{
		(*it)->solveTNF(_grid.getGrid(), dt);
	}
}


void GrSimulation::updateIL10Dynamics(double dt)
{
	for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
	{
		(*it)->solveIL10(_grid.getGrid(), dt);
	}
	for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
	{
		(*it)->solveIL10(_grid.getGrid(), dt);
	}
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
		(*it)->solveIL10(_grid.getGrid(), dt);
	}
	for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
	{
		(*it)->solveIL10(_grid.getGrid(), dt);
	}
}



void GrSimulation::updateTNFandIL10Dynamics(double dt)
{
	for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
	{
		(*it)->solveTNFandIL10(_grid.getGrid(), dt);
	}
	for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
	{
		(*it)->solveTNFandIL10(_grid.getGrid(), dt);
	}
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
		(*it)->solveTNFandIL10(_grid.getGrid(), dt);
	}
	for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
	{
		(*it)->solveTNFandIL10(_grid.getGrid(), dt);
	}
}


void GrSimulation::updateNFkBandTNFDynamics(double dt)
{
	for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
	{
		for (int i = 0; i < _PARAM(PARAM_GR_NF_KB_TIME_COEFF); i++)
		{
			(*it)->solveNFkBandTNF(_grid.getGrid(), dt/_PARAM(PARAM_GR_NF_KB_TIME_COEFF));
		}
	}
	for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
	{
		(*it)->solveTNF(_grid.getGrid(), dt);
	}
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
		(*it)->solveTNF(_grid.getGrid(), dt);
	}
	for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
	{
		(*it)->solveTNF(_grid.getGrid(), dt);
	}
}


void GrSimulation::updateNFkBandTNFandIL10Dynamics(double dt)
{
	for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
	{
		for (int i = 0; i < _PARAM(PARAM_GR_NF_KB_TIME_COEFF); i++)
		{
			(*it)->solveTNFandIL10andNFkB(_grid.getGrid(), dt/_PARAM(PARAM_GR_NF_KB_TIME_COEFF));
		}
	}
	for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
	{
		(*it)->solveTNFandIL10(_grid.getGrid(), dt);
	}
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
		(*it)->solveTNFandIL10(_grid.getGrid(), dt);
	}
	for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
	{
		(*it)->solveTNFandIL10(_grid.getGrid(), dt);
	}
}
#endif


void GrSimulation::moveMacrophages()
{
	const int timeResting = _PARAM(PARAM_MAC_MOVEMENT_RESTING);
	const int timeInfected = _PARAM(PARAM_MAC_MOVEMENT_INFECTED);
	const int timeActive = _PARAM(PARAM_MAC_MOVEMENT_ACTIVE);
	
	for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
	{
		// move resting macrophages
		if (_time % timeResting == 0 && (*it)->getState() == Mac::MAC_RESTING)
		{
			(*it)->move(_grid.getGrid());
		}
		// move active macrophages
		else if (_time % timeActive == 0 && (*it)->getState() == Mac::MAC_ACTIVE)
		{
			(*it)->move(_grid.getGrid());
		}
		// move infected macrophages
		else if (_time % timeInfected == 0 && (*it)->getState() == Mac::MAC_INFECTED)
		{
			(*it)->move(_grid.getGrid());
		}
	} 
}

void GrSimulation::moveTcells()
{
		
	for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
	{
		(*it)->move(_grid.getGrid());
	}
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
		(*it)->move(_grid.getGrid());
	}
	for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
	{
		(*it)->move(_grid.getGrid());
	} 
}

void GrSimulation::growExtMtb()
{
	const double growthRate = _PARAM(PARAM_EXTMTB_GROWTH_RATE) - 1;
	const double upperBound = _PARAM(PARAM_EXTMTB_UPPER_BOUND);

  GrGrid& g = _grid.getGrid();
  const Pos& dim = g.getRange();
  Pos p;
	for (p.x = 0; p.x < dim.x; p.x++)
	{
		for (p.y = 0; p.y < dim.y; p.y++)
		{
			Scalar& extMtb = g.extMTB(p);
			
			if (g.isCaseated(p))
			{
				// Bacteria don't gp.x in caseated compartments
				++_stats.getNrCaseated();
				_stats.getTotNonRepExtMtb()+=(extMtb);
				_stats.getTotExtMtb()+=(extMtb);
			}
			else
			{
				double dExtMtb = growthRate * extMtb * (1.0 - extMtb / upperBound);
        extMtb += dExtMtb;
				_stats.getTotExtMtb() += (extMtb);
			}

			_stats.getTotMacAttractant() += (g.macAttractant(p));
			_stats.getTotTNF() += (g.TNF(p));
      _stats.getTotIL10() += (g.il10(p));
			_stats.getTotCCL2() += (g.CCL2(p));
			_stats.getTotCCL5() += (g.CCL5(p));
			_stats.getTotCXCL9() += (g.CXCL9(p));

			if (g.TNF(p) >= _areaThreshold)
				++_stats.getAreaTNF();

			if (g.getCellDensity(p) >= _areaThresholdCellDensity)
				++_stats.getAreaCellDensity();
		}
	}
}

void GrSimulation::adjustTNFDegradation(double dt)
{	
  Pos pos;
  GrGrid& grid = _grid.getGrid();
	for (pos.x = 0; pos.x < grid.getRange().x; pos.x++)
	{
		for (pos.y = 0; pos.y < grid.getRange().y; pos.y++)
		{
			// simulate the effect of TNF internalization by cells in the form of degradation. only for TNF
			double dtnf;
			Scalar tnf = grid.TNF(pos);
			if (grid.hasAgentType(MAC, pos))
			{
				dtnf = -_PARAM(PARAM_GR_K_INT1) * (tnf / (tnf + _PARAM(PARAM_GR_KD1) * 48.16e11)) * 1500 * dt * 0.4;
				tnf += dtnf;
			}
            
            if (grid.hasTcell(pos))
            {
                dtnf = -_PARAM(PARAM_GR_K_INT1) * (tnf / (tnf + _PARAM(PARAM_GR_KD1) * 48.16e11)) * 800 * dt * 0.4;
                tnf += dtnf;  
            }
            
			// dtnf = -_PARAM(PARAM_GR_K_INT1) * (tnf / (tnf + _PARAM(PARAM_GR_KD1) * 48.16e11)) * 800 * (cell.hasTcell()) * dt * 0.4;
			// tnf += dtnf;
		}
	}
}

void GrSimulation::adjustFauxDegradation(double dt)
{
    
  for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
  {
      (*it)->solveDegradation(_grid.getGrid(), dt, _tnfrDynamics, _il10rDynamics);
  }
  for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
	{
		(*it)->solveDegradation(_grid.getGrid(), dt, _tnfrDynamics, _il10rDynamics);
	}
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
		(*it)->solveDegradation(_grid.getGrid(), dt, _tnfrDynamics, _il10rDynamics);
	}
	for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
	{
		(*it)->solveDegradation(_grid.getGrid(), dt, _tnfrDynamics, _il10rDynamics);
	}
}

void GrSimulation::updateT_Test()
{
	for (int i = 0; i < NOUTCOMES; i++)
	{
		if (_pTTest[i])
		{
			_pTTest[i]->update(_time, i, _stats);
		}
	}
}

void GrSimulation::performT_Test()
{
	for (int i = 0; i < NOUTCOMES; i++)
	{
		if (_pTTest[i])
		{
			_pTTest[i]->perform(i, _stats);
		}
	}
}

void GrSimulation::setDiffusionMethod(DiffusionMethod method)
{
	if (_pDiffusion->getMethod() != method)
	{
		delete _pDiffusion;

		switch (method)
		{
		case DIFF_SOR_CORRECT:
			_pDiffusion = new GrDiffusionBTCS();
			break;
		case DIFF_SOR_WRONG:
			_pDiffusion = new GrDiffusionWrongBTCS();
			break;
    case DIFF_REC_EQ:
        std::clog<<"Warning: DIFF_REC_EQ disabled, defaulting to DIFF_REC_EQ_SWAP"<<std::endl;
		case DIFF_REC_EQ_SWAP:
			_pDiffusion = new GrDiffusionFTCS_Swap();
			break;
    case DIFF_ADE_SWAP:
        _pDiffusion = new GrDiffusionADE_Swap();
        break;
                
		}
    if(method != DIFF_ADE_SWAP && _PARAM(PARAM_GR_DT_DIFFUSION) > 12) {
     std::clog<<("*** WARNING: This diffusion method is unstable for timesteps greater than 12 seconds")<<std::endl;
    }
	}
}

RecruitmentMethod GrSimulation::getRecruitmentMethod()
{
	return _pRecruitment->getMethod();
}

void GrSimulation::setRecruitmentMethod(RecruitmentMethod method)
{
	if (_pRecruitment)
	{
		if (method == _pRecruitment->getMethod())
		{
			// Already have the correct type of recruitment object.
			// Don't want to replace it since it may have been deserialized from a saved simulation state.
			return;
		}

		delete _pRecruitment;
	}

  switch (method)
  {
  case RECR_PROB:
	  _pRecruitment = new RecruitmentProb();
    break;

  case RECR_LN_ODE_PURE:
	  _pRecruitment = new RecruitmentLnODEPure();
      break;

  default:
    std::cerr << "Invalid recruitment method: " << method << std::endl;
    exit(1);
  }
}

void GrSimulation::deserializeRecruitmentMethod(RecruitmentMethod method, std::istream& in)
{
  assert(in.good());
	if (_pRecruitment)
		delete _pRecruitment;

  switch (method)
  {
  case RECR_PROB:
	  _pRecruitment = new RecruitmentProb();
    break;

  case RECR_LN_ODE_PURE:
	  _pRecruitment = new RecruitmentLnODEPure(in);
      break;

  default:
	    std::cerr << "Invalid recruitment method: " << method << std::endl;
	    exit(1);
  }
  assert(in.good());
}

void GrSimulation::setOutcomeMethod(int index, OutcomeMethod method, double alpha,
	int testPeriod, int samplePeriod)
{
	assert(0 <= index && index < NOUTCOMES);

	if (_pTTest[index] && _pTTest[index]->getMethod() == method)
	{
		_pTTest[index]->setAlpha(alpha);
		_pTTest[index]->setTestPeriod(testPeriod);
		_pTTest[index]->setSamplePeriod(samplePeriod);
	}
	else
	{
		delete _pTTest[index];

		switch (method)
		{
		case OUTCOME_AREA:
			_pTTest[index] = new AreaTest(alpha, testPeriod, samplePeriod);
			break;
		case OUTCOME_MTB:
			_pTTest[index] = new MtbTest(alpha, testPeriod, samplePeriod);
			break;
		case OUTCOME_NONE:
			_pTTest[index] = NULL;
			break;
		}
	}
}

void GrSimulation::shuffleCells()
{
#if 1
  std::random_shuffle(_macList.begin(), _macList.end(), g_Rand);
  std::random_shuffle(_tgamList.begin(), _tgamList.end(), g_Rand);
  std::random_shuffle(_tcytList.begin(), _tcytList.end(), g_Rand);
  std::random_shuffle(_tregList.begin(), _tregList.end(), g_Rand);
#else
  //TODO: replace with random_shuffle
	MacList::iterator mac1 = _macList.begin();
	MacList::iterator mac2 = _macList.end();
	TgamList::iterator tgam1 = _tgamList.begin();
	TgamList::iterator tgam2 = _tgamList.end();
	TcytList::iterator tcyt1 = _tcytList.begin();
	TcytList::iterator tcyt2 = _tcytList.end();
	TregList::iterator treg1 = _tregList.begin();
	TregList::iterator treg2 = _tregList.end();
	MacList::size_type m = 1;
	MacList::size_type n = 1;
	
	// Don't shuffle an empty list.
	// The list size is unsigned so the shuffle loop limit, _list size - 1, is also unsigned.
	// When list size is 0, the limit is a large positive number instead of -1.

	if (_macList.size() > 0)
	{
		while (m + n < _macList.size() - 1)
		{
			if (g_Rand.getReal() < 0.5)
			{
				mac1++;
				m++;
			}
			if (g_Rand.getReal() < 0.5)
			{
				mac2--;
				n++;
			}
			std::swap(mac1, mac2);
		}
	}
	if (_tgamList.size() > 0)
	{
		m = 1;
		n = 1;
		while (m + n < _tgamList.size() - 1)
		{
			if (g_Rand.getReal() < 0.5)
			{
				tgam1++;
				m++;
			}
			if (g_Rand.getReal() < 0.5)
			{
				tgam2--;
				n++;
			}
			std::swap(tgam1, tgam2);
		}	
	}
	
	if (_tcytList.size() > 0)
	{
		m = 1;
		n = 1;
		while (m + n < _tcytList.size() - 1)
		{
			if (g_Rand.getReal() < 0.5)
			{
				tcyt1++;
				m++;
			}
			if (g_Rand.getReal() < 0.5)
			{
				tcyt2--;
				n++;
			}
			std::swap(tcyt1, tcyt2);
		}
	}
	
	if (_tregList.size() > 0)
	{
		m = 1;
		n = 1;
		while (m + n < _tregList.size() - 1)
		{
			if (g_Rand.getReal() < 0.5)
			{
				treg1++;
				m++;
			}
			if (g_Rand.getReal() < 0.5)
			{
				treg2--;
				n++;
			}
			std::swap(treg1, treg2);
		}
	}	
#endif
}
