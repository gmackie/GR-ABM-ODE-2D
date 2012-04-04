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
#include "grdiffusionftcsswap.h"
#include "areatest.h"
#include "mtbtest.h"
#include "recruitmentprob.h"
#include "serialization.h"

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
	, _pDiffusion(new GrDiffusionWrongBTCS())
	, _pTTest()
	, _pRecruitment(NULL)
	, _tnfrDynamics(false)
	, _nfkbDynamics(false)
    , _il10rDynamics(false)
    , _tgammatransition(false)
	, _tnfDepletionTimeStep(-1)
    , _il10DepletionTimeStep(-1)
	, _tcellRecruitmentBegun(false)
{
	for (int i = 0; i < NOUTCOMES; i++)
		_pTTest[i] = NULL;
}

GrSimulation::~GrSimulation()
{
	delete _pDiffusion;
	for (int i = 0; i < NOUTCOMES; i++)
		delete _pTTest[i];
}

void GrSimulation::serialize(std::ostream& out) const
{
	assert(out.good());

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

	out << _tnfrDynamics << std::endl;
	out << _nfkbDynamics << std::endl;
    out << _il10rDynamics << std::endl;
    out << _tgammatransition << std::endl;
	out << _tnfDepletionTimeStep << std::endl;
    out << _il10DepletionTimeStep << std::endl;
    out << _tcellRecruitmentBegun <<std::endl;

	// serialize grid
	_grid.serialize(out);

	// Serialize the agent class.
	Agent::classSerialize(out);

	// serialize macs
	out << _macList.size() << std::endl;
	for (MacList::const_iterator it = _macList.begin(); it != _macList.end(); it++)
	{
		it->serialize(out);
	}

	// serialize tgam cells
	out << _tgamList.size() << std::endl;
	for (TgamList::const_iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
	{
		it->serialize(out);
	}

	// serialize tcyt cells
	out << _tcytList.size() << std::endl;
	for (TcytList::const_iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
		it->serialize(out);
	}

	// serialize treg cells
	out << _tregList.size() << std::endl;
	for (TregList::const_iterator it = _tregList.begin(); it != _tregList.end(); it++)
	{
		it->serialize(out);
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

	in >> _tnfrDynamics;
	in >> _nfkbDynamics;
    in >> _il10rDynamics;
    in >> _tgammatransition;
	in >> _tnfDepletionTimeStep;
    in >> _il10DepletionTimeStep;
    in >> _tcellRecruitmentBegun;

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
		// create a dummy mac
		_macList.push_back(Mac());
		Mac* pMac = &_macList.back();

		// update attributes of dummy mac and add to grid
		pMac->deserialize(in);
		bool addedAgent = _grid.getGrid().addAgent(pMac, pMac->getPosition());
		if (!addedAgent)
		{
			cerr << "Error deserializing mac " << i << " at grid position " << pMac->getPosition() << ":addAgent failed."<< endl;
			exit(1);
		}
	}

	// deserialize tgam cells
	_tgamList.clear();
	in >> intVal;
	for (int i = 0; i < intVal; i++)
	{
		// create a dummy tgam cell
		_tgamList.push_back(Tgam());
		Tgam* pTgam = &_tgamList.back();

		// update attributes of dummy tgam and add to grid
		pTgam->deserialize(in);
		bool addedAgent = _grid.getGrid().addAgent(pTgam, pTgam->getPosition());
		if (!addedAgent)
		{
			cerr << "Error deserializing tgam " << i << " at grid position " << pTgam->getPosition() << ":addAgent failed."<< endl;
			exit(1);
		}

	}

	// deserialize tcyt cells
	_tcytList.clear();
	in >> intVal;
	for (int i = 0; i < intVal; i++)
	{
		// create a dummy tcyt cell
		_tcytList.push_back(Tcyt());
		Tcyt* pTcyt = &_tcytList.back();

		// update attributes of dummy tcyt and add to grid
		pTcyt->deserialize(in);
		bool addedAgent = _grid.getGrid().addAgent(pTcyt, pTcyt->getPosition());
		if (!addedAgent)
		{
			cerr << "Error deserializing tcyt " << i << " at grid position " << pTcyt->getPosition() << ":addAgent failed."<< endl;
			exit(1);
		}

	}

	// deserialize treg cells
	_tregList.clear();
	in >> intVal;
	for (int i = 0; i < intVal; i++)
	{
		// create a dummy treg cell
		_tregList.push_back(Treg());
		Treg* pTreg = &_tregList.back();

		// update attributes of dummy mac and add to grid
		pTreg->deserialize(in);
		bool addedAgent = _grid.getGrid().addAgent(pTreg, pTreg->getPosition());
		if (!addedAgent)
		{
			cerr << "Error deserializing treg " << i << " at grid position " << pTreg->getPosition() << ":addAgent failed."<< endl;
			exit(1);
		}

	}

	// deserialize statistics
	_statsPrevious.deserialize(in);
	_stats.deserialize(in);

	// deserialize random number generator
	g_Rand.deserialize(in);

	if (!Serialization::readFooter(in, GrSimulation::_ClassName))
	{
		exit(1);
	}
}

void GrSimulation::init(Scalar molecularTrackingRadius)
{
	// initialize the sources
	_grid.initSources();

	PosVector initMacs = Params::getInstance()->getInitialMacs();
	PosVector initExtMtb = Params::getInstance()->getInitialExtMtb();
	const int maxMacAge = _PARAM(PARAM_MAC_AGE);

	// Place initial infected macrophages on the grid
	for (PosVector::const_iterator it = initMacs.begin(); it != initMacs.end(); it++)
	{
		Mac* pMac = createMac(it->x, it->y, g_Rand.getInt(-1, -maxMacAge), MAC_INFECTED, false, false);
		
		if (_nfkbDynamics)
		{
			// initialize NF-kB signaling from steady-state (12 hours equilibrium)
			for (int i = 0; i < 7200*_PARAM(PARAM_GR_NF_KB_TIME_COEFF); ++i)
				pMac->solveNFkBODEsEquilibrium(6.00/_PARAM(PARAM_GR_NF_KB_TIME_COEFF));
		}

		pMac->setIntMtb(1);
		_statsPrevious.incTotIntMtb(1);
		_stats.incTotIntMtb(1);
	}

	// Place initial extracellular bacteria on the grid
	for (PosVector::const_iterator it = initExtMtb.begin(); it != initExtMtb.end(); it++)
	{
		_grid.getGrid().extMTB(*it) += 1;
		_statsPrevious.incTotExtMtb(1);
		_stats.incTotExtMtb(1);
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
			Mac* newMac = createMac(row, col, g_Rand.getInt(-1, -maxMacAge), MAC_RESTING, false, false);
			if (_nfkbDynamics)
			{
				// initialize NF-kB signaling from steady-state (12 hours equilibrium)
				for (int i = 0; i < 7200*_PARAM(PARAM_GR_NF_KB_TIME_COEFF); ++i)
					newMac->solveNFkBODEsEquilibrium(6.00/_PARAM(PARAM_GR_NF_KB_TIME_COEFF));
			}
			count--;
		}
	}

	initMolecularTracking(molecularTrackingRadius);
}

void GrSimulation::initMolecularTracking(Scalar molecularTrackingRadius)
{
	const Pos center = _grid.getCenter();

	for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
	{
		if (it->getPosition().distance(center) <= molecularTrackingRadius)
		{
			it->setTrackMolecularDynamics(true);
		}
	}

	for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
	{
		if (it->getPosition().distance(center) <= molecularTrackingRadius)
		{
			it->setTrackMolecularDynamics(true);
		}
	}

	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
		if (it->getPosition().distance(center) <= molecularTrackingRadius)
		{
			it->setTrackMolecularDynamics(true);
		}
	}

	for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
	{
		if (it->getPosition().distance(center) <= molecularTrackingRadius)
		{
			it->setTrackMolecularDynamics(true);
		}
	}
}

void GrSimulation::solve()
{
	// here we perform a timestep, which is 10 minutes
	_time++;

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
    
	// calculate diffusion every 10 minutes
	// dt = 6s, solve for 10 minutes = 100 * 6 seconds
	double dt = 6;
	for (int t = 0; t < 100; t++) 
	{
		secreteFromMacrophages(tnfDepletion, il10Depletion);
		secreteFromTcells(tnfDepletion, il10Depletion);
		secreteFromCaseations();
		_pDiffusion->diffuse(_grid);

        if (_nfkbDynamics)
		{
			if (_il10rDynamics) {
                updateNFkBandTNFandIL10Dynamics(dt);
            }
			else
			{
				updateNFkBandTNFDynamics(dt);
			}
		}
		else if (_tnfrDynamics && _il10rDynamics)
        {
            updateTNFandIL10Dynamics(dt, t);
        }
        
        else if (_tnfrDynamics && !_il10rDynamics)
		{
			updateTNFDynamics(dt);
		}
		
        else if (_il10rDynamics && !_tnfrDynamics)
        {
            updateIL10Dynamics(dt);
        }
        
        else if (!_il10rDynamics || !_tnfrDynamics)
		{
			adjustFauxDegradation(dt);
            //adjustTNFDegradation(dt);
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
}

void GrSimulation::updateStates()
{
	// update states and remove dead macrophages from lists
	for (MacList::iterator it = _macList.begin(); it != _macList.end();)
	{
		Mac& mac = *it;
		mac.updateState();
		_stats.updateAgentStatistics(&mac);

		if (!mac.isDead())
			_stats.incTotIntMtb(mac.getIntMtb());

		if (mac.getNFkB())
			_stats.updateMacNFkBStatistics((MacState)mac.getState());

		if (mac.getStat1())
			_stats.updateMacStat1Statistics((MacState)mac.getState());

		if (mac.isDeactivated())
			_stats.updateMacDeactStatistics((MacState)mac.getState());

		if (mac.isDead())
		{
			assert_res(_grid.getGrid().removeAgent(&mac));
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
		Tgam& tgam = *it;
		tgam.updateState();
		_stats.updateAgentStatistics(&tgam);

		if (tgam.isDead())
		{
			assert_res(_grid.getGrid().removeAgent(&tgam));
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
		Tcyt& tcyt = *it;
		tcyt.updateState();
		_stats.updateAgentStatistics(&tcyt);

		if (tcyt.isDead())
		{
			assert_res(_grid.getGrid().removeAgent(&tcyt));
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
		Treg& treg = *it;
		treg.updateState();
		_stats.updateAgentStatistics(&treg);

		if (it->isDead())
		{
			assert_res(_grid.getGrid().removeAgent(&treg));
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
		it->computeNextState(_time, _grid.getGrid(), _stats, _tnfrDynamics, _nfkbDynamics, _il10rDynamics, _tgammatransition);
	}
	for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
	{
		it->computeNextState(_time, _grid.getGrid(), _stats, _tnfrDynamics, _nfkbDynamics, _il10rDynamics, _tgammatransition);
	}
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
		it->computeNextState(_time, _grid.getGrid(), _stats, _tnfrDynamics, _nfkbDynamics, _il10rDynamics, _tgammatransition);
	}
	for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
	{
		it->computeNextState(_time, _grid.getGrid(), _stats, _tnfrDynamics, _nfkbDynamics, _il10rDynamics, _tgammatransition);
	} 
}

void GrSimulation::secreteFromMacrophages(bool tnfDepletion, bool il10Depletion)
{
	for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
	{
		it->secrete(_grid.getGrid(), _tnfrDynamics, _nfkbDynamics, tnfDepletion, _il10rDynamics, il10Depletion);
	}
}

void GrSimulation::secreteFromTcells(bool tnfDepletion, bool il10Depletion)
{
	for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
	{
		it->secrete(_grid.getGrid(), _tnfrDynamics, _nfkbDynamics, tnfDepletion, _il10rDynamics, il10Depletion);
	}
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
		it->secrete(_grid.getGrid(), _tnfrDynamics, _nfkbDynamics, tnfDepletion, _il10rDynamics, il10Depletion);
	}
    for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
	{
		it->secrete(_grid.getGrid(), _tnfrDynamics, _nfkbDynamics, tnfDepletion, _il10rDynamics, il10Depletion);
	}
}

void GrSimulation::secreteFromCaseations()
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
				g.incCCL2(p, (0.25 * _PARAM(PARAM_MAC_SEC_RATE_CCL2)));
				g.incCCL5(p, (0.25 * _PARAM(PARAM_MAC_SEC_RATE_CCL5)));
				g.incCXCL9(p,  (0.25 * _PARAM(PARAM_MAC_SEC_RATE_CXCL9)));
				g.incmacAttractant(p, (_PARAM(PARAM_GR_SEC_RATE_ATTRACTANT)));
			}
		}
	}
}

void GrSimulation::updateTNFDynamics(double dt)
{
	for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
	{
		it->solveTNF(_grid.getGrid(), dt);
	}
	for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
	{
		it->solveTNF(_grid.getGrid(), dt);
	}
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
		it->solveTNF(_grid.getGrid(), dt);
	}
	for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
	{
		it->solveTNF(_grid.getGrid(), dt);
	}
}


void GrSimulation::updateIL10Dynamics(double dt)
{
	for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
	{
		it->solveIL10(_grid.getGrid(), dt);
	}
	for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
	{
		it->solveIL10(_grid.getGrid(), dt);
	}
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
		it->solveIL10(_grid.getGrid(), dt);
	}
	for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
	{
		it->solveIL10(_grid.getGrid(), dt);
	}
}



void GrSimulation::updateTNFandIL10Dynamics(double dt, double currenttime)
{
	for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
	{
		it->solveTNFandIL10(_grid.getGrid(), _stats, dt, currenttime);
	}
	for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
	{
		it->solveTNFandIL10(_grid.getGrid(), _stats, dt, currenttime);
	}
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
		it->solveTNFandIL10(_grid.getGrid(), _stats, dt, currenttime);
	}
	for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
	{
		it->solveTNFandIL10(_grid.getGrid(), _stats, dt, currenttime);
	}
}


void GrSimulation::updateNFkBandTNFDynamics(double dt)
{
	for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
	{
		for (int i = 0; i < _PARAM(PARAM_GR_NF_KB_TIME_COEFF); i++)
		{
			it->solveNFkBandTNF(_grid.getGrid(), dt/_PARAM(PARAM_GR_NF_KB_TIME_COEFF));
		}
	}
	for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
	{
		it->solveTNF(_grid.getGrid(), dt);
	}
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
		it->solveTNF(_grid.getGrid(), dt);
	}
	for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
	{
		it->solveTNF(_grid.getGrid(), dt);
	}
}


void GrSimulation::updateNFkBandTNFandIL10Dynamics(double dt)
{
	for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
	{
		for (int i = 0; i < _PARAM(PARAM_GR_NF_KB_TIME_COEFF); i++)
		{
			it->solveTNFandIL10andNFkB(_grid.getGrid(), dt/_PARAM(PARAM_GR_NF_KB_TIME_COEFF));
		}
	}
	for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
	{
		it->solveTNFandIL10(_grid.getGrid(), _stats, dt, 0.0);
	}
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
		it->solveTNFandIL10(_grid.getGrid(), _stats, dt, 0.0);
	}
	for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
	{
		it->solveTNFandIL10(_grid.getGrid(), _stats, dt, 0.0);
	}
}



void GrSimulation::moveMacrophages()
{
	const int timeResting = _PARAM(PARAM_MAC_MOVEMENT_RESTING);
	const int timeInfected = _PARAM(PARAM_MAC_MOVEMENT_INFECTED);
	const int timeActive = _PARAM(PARAM_MAC_MOVEMENT_ACTIVE);
	
	for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
	{
		// move resting macrophages
		if (_time % timeResting == 0 && it->getState() == MAC_RESTING)
		{
			it->move(_grid.getGrid());
		}
		// move active macrophages
		else if (_time % timeActive == 0 && it->getState() == MAC_ACTIVE)
		{
			it->move(_grid.getGrid());
		}
		// move infected macrophages
		else if (_time % timeInfected == 0 && it->getState() == MAC_INFECTED)
		{
			it->move(_grid.getGrid());
		}
	} 
}

void GrSimulation::moveTcells()
{
		
	for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
	{
		it->move(_grid.getGrid());
	}
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
		it->move(_grid.getGrid());
	}
	for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
	{
		it->move(_grid.getGrid());
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
				_stats.incNrCaseated();
				_stats.incTotNonRepExtMtb(extMtb);
				_stats.incTotExtMtb(extMtb);
			}
			else
			{
				double dExtMtb = growthRate * extMtb * (1.0 - extMtb / upperBound);
        extMtb += dExtMtb;
				_stats.incTotExtMtb(extMtb);
			}

			_stats.incTotMacAttractant(g.macAttractant(p));
			_stats.incTotTNF(g.TNF(p));
            _stats.incTotIL10(g.il10(p));
			_stats.incTotCCL2(g.CCL2(p));
			_stats.incTotCCL5(g.CCL5(p));
			_stats.incTotCXCL9(g.CXCL9(p));

			if (g.TNF(p) >= _areaThreshold)
				_stats.incAreaTNF();

			if (g.getCellDensity(p) >= _areaThresholdCellDensity)
			{
				_stats.incAreaCellDensity();
			}
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
        it->solveDegradation(_grid.getGrid(), dt, _tnfrDynamics, _il10rDynamics);
    }
    for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
	{
		it->solveDegradation(_grid.getGrid(), dt, _tnfrDynamics, _il10rDynamics);
	}
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
		it->solveDegradation(_grid.getGrid(), dt, _tnfrDynamics, _il10rDynamics);
	}
	for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
	{
		it->solveDegradation(_grid.getGrid(), dt, _tnfrDynamics, _il10rDynamics);
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
		}
	}
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
}
