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
#include "grdiffusionftcs.h"
#include "grdiffusionftcsswap.h"
#include "areatest.h"
#include "mtbtest.h"
#include "recruitmentprob.h"
#include "serialization.h"

const std::string GrSimulation::_ClassName = "GrSimulation";

GrSimulation::GrSimulation()
	: _time(0)
	, _grid()
	, _macList()
	, _tgamList()
	, _tcytList()
	, _tregList()
	, _stats()
	, _areaThreshold(0.5f)
	, _areaThresholdCellDensity(0.5f)
	, _pDiffusion(new GrDiffusionWrongBTCS())
	, _pTTest()
	, _pRecruitment(NULL)
	, _tnfrDynamics(false)
	, _nfkbDynamics(false)
    , _il10rDynamics(false)
	, _tnfDepletionTimeStep(-1)
    , _il10DepletionTimeStep(-1)
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
	out << _tnfDepletionTimeStep << std::endl;
    out << _il10DepletionTimeStep << std::endl;

	// serialize grid
	_grid.serialize(out);

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
	in >> _tnfDepletionTimeStep;
    in >> _il10DepletionTimeStep;

	// deserialize grid
	_grid.deserialize(in);

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
		assert_res(_grid(pMac->getRow(), pMac->getCol()).addAgent(pMac));
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
		assert_res(_grid(pTgam->getRow(), pTgam->getCol()).addAgent(pTgam));
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
		assert_res(_grid(pTcyt->getRow(), pTcyt->getCol()).addAgent(pTcyt));
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
		assert_res(_grid(pTreg->getRow(), pTreg->getCol()).addAgent(pTreg));
	}

	// deserialize statistics
	_stats.deserialize(in);

	// deserialize random number generator
	g_Rand.deserialize(in);

	if (!Serialization::readFooter(in, GrSimulation::_ClassName))
	{
		exit(1);
	}
}

void GrSimulation::init()
{
	// initialize the sources
	_grid.initSources();

	PosVector initMacs = Params::getInstance()->getInitialMacs();
	PosVector initExtMtb = Params::getInstance()->getInitialExtMtb();
	const int maxMacAge = _PARAM(PARAM_MAC_AGE);

	// Place initial infected macrophages on the grid
	for (PosVector::const_iterator it = initMacs.begin(); it != initMacs.end(); it++)
	{
		Mac* pMac = createMac(it->first, it->second,
			g_Rand.getInt(-1, -maxMacAge), MAC_INFECTED, false, false);
		
		if (_nfkbDynamics)
		{
			// initialize NF-kB signaling from steady-state (12 hours equilibrium)
			for (int i = 0; i < 7200*_PARAM(PARAM_GR_NF_KB_TIME_COEFF); ++i)
				pMac->solveNFkBODEsEquilibrium(6.00/_PARAM(PARAM_GR_NF_KB_TIME_COEFF));
		}

		pMac->setIntMtb(1);
		_stats.incTotIntMtb(1);
	}

	// Place initial extracellular bacteria on the grid
	for (PosVector::const_iterator it = initExtMtb.begin(); it != initExtMtb.end(); it++)
	{
		_grid(it->first, it->second).incExtMtb(1);
		_stats.incTotExtMtb(1);
	}

	// Add a fixed number (PARAM_MAC_INIT_NUMBER) of resting
	// macrophages to the grid
	int count = _PARAM(PARAM_MAC_INIT_NUMBER);

	int nrGridCompartments = (NROWS*NCOLS);
	if (count > nrGridCompartments)
	{
		std::cerr << "The number of initial resting macrophages to place on the grid, " << count << ", is > the number of grid compartments, " << nrGridCompartments << "." << std::endl;
		exit(1);
	}

	while (count > 0)
	{
		int row = g_Rand.getInt(NROWS);
		int col = g_Rand.getInt(NCOLS);
		if (!_grid(row, col).hasMac())
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
            updateNFkBandTNFDynamics(dt);
		}
		else if (_tnfrDynamics && _il10rDynamics)
        {
            updateTNFandIL10Dynamics(dt);
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
			adjustFauxDegradation(dt, _tnfrDynamics, _il10rDynamics);
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
	
	// shuffle order of cells every hour for randomization in movement
	if (_time % 6 == 0)
	{
		shuffleCells();
	}

	// perform t-test
	updateT_Test();
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
			_stats.updateMacNFkBStatistics(mac.getState());

		if (mac.getStat1())
			_stats.updateMacStat1Statistics(mac.getState());

		if (mac.isDeactivated())
			_stats.updateMacDeactStatistics(mac.getState());

		if (mac.isDead())
		{
			assert_res(_grid(mac.getRow(), mac.getCol()).removeAgent(&mac));
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
			assert_res(_grid(tgam.getRow(), tgam.getCol()).removeAgent(&tgam));
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
			assert_res(_grid(tcyt.getRow(), tcyt.getCol()).removeAgent(&tcyt));
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
			assert_res(_grid(treg.getRow(), treg.getCol()).removeAgent(&treg));
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
		it->computeNextState(_time, _grid.getGrid(), _stats, _tnfrDynamics, _nfkbDynamics, _il10rDynamics);
	}
	for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
	{
		it->computeNextState(_time, _grid.getGrid(), _stats, _tnfrDynamics, _nfkbDynamics, _il10rDynamics);
	}
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
		it->computeNextState(_time, _grid.getGrid(), _stats, _tnfrDynamics, _nfkbDynamics, _il10rDynamics);
	}
	for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
	{
		it->computeNextState(_time, _grid.getGrid(), _stats, _tnfrDynamics, _nfkbDynamics, _il10rDynamics);
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
	for (int i = 0; i < NROWS; i++)
	{
		for (int j = 0; j < NCOLS; j++)
		{
			GridCell& cell = _grid(i, j);
			if (cell.isCaseated())
			{
				cell.incCCL2(0.25 * _PARAM(PARAM_MAC_SEC_RATE_CCL2));
				cell.incCCL5(0.25 * _PARAM(PARAM_MAC_SEC_RATE_CCL5));
				cell.incCXCL9(0.25 * _PARAM(PARAM_MAC_SEC_RATE_CXCL9));
			    cell.incMacAttractant(_PARAM(PARAM_GR_SEC_RATE_ATTRACTANT));
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



void GrSimulation::updateTNFandIL10Dynamics(double dt)
{
	for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
	{
		it->solveTNFandIL10(_grid.getGrid(), dt);
	}
	for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
	{
		it->solveTNFandIL10(_grid.getGrid(), dt);
	}
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
		it->solveTNFandIL10(_grid.getGrid(), dt);
	}
	for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
	{
		it->solveTNFandIL10(_grid.getGrid(), dt);
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
		it->solveTNFandIL10(_grid.getGrid(), dt);
	}
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
		it->solveTNFandIL10(_grid.getGrid(), dt);
	}
	for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
	{
		it->solveTNFandIL10(_grid.getGrid(), dt);
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

	for (int row = 0; row < NROWS; row++)
	{
		for (int col = 0; col < NCOLS; col++)
		{
			GridCell& cell = _grid(row, col);
			
			double extMtb = cell.getExtMtb();
			
			if (cell.isCaseated())
			{
				// Bacteria don't grow in caseated compartments
				_stats.incNrCaseated();
				_stats.incTotNonRepExtMtb(extMtb);
				_stats.incTotExtMtb(extMtb);
			}
			else
			{
				double dExtMtb = growthRate * extMtb * (1.0 - extMtb / upperBound);
				cell.incExtMtb(dExtMtb);
				_stats.incTotExtMtb(cell.getExtMtb());
			}

			_stats.incTotMacAttractant(cell.getMacAttractant());
			_stats.incTotTNF(cell.getTNF());
            _stats.incTotIL10(cell.getIL10());
			_stats.incTotCCL2(cell.getCCL2());
			_stats.incTotCCL5(cell.getCCL5());
			_stats.incTotCXCL9(cell.getCXCL9());

			if (cell.getTNF() >= _areaThreshold)
				_stats.incAreaTNF();

			if (cell.getCellDensity(_grid.getGrid()) >= _areaThresholdCellDensity)
			{
				_stats.incAreaCellDensity();
			}
		}
	}
}

void GrSimulation::adjustTNFDegradation(double dt)
{	
	for (int row = 0; row < NROWS; row++)
	{
		for (int col = 0; col < NCOLS; col++)
		{
			GridCell& cell = _grid(row, col);
			
			// simulate the effect of TNF internalization by cells in the form of degradation. only for TNF
			double dtnf;
			double tnf = cell.getTNF();
			if (cell.hasMac())
			{
				dtnf = -_PARAM(PARAM_GR_K_INT1) * (tnf / (tnf + _PARAM(PARAM_GR_KD1) * 48.16e11)) * 1500 * dt * 0.4;
				tnf += dtnf;
			}
            
            if (cell.hasTcell())
            {
                dtnf = -_PARAM(PARAM_GR_K_INT1) * (tnf / (tnf + _PARAM(PARAM_GR_KD1) * 48.16e11)) * 800 * dt * 0.4;
                tnf += dtnf;  
            }
            
			// dtnf = -_PARAM(PARAM_GR_K_INT1) * (tnf / (tnf + _PARAM(PARAM_GR_KD1) * 48.16e11)) * 800 * (cell.hasTcell()) * dt * 0.4;
			// tnf += dtnf;
			
			cell.setTNF(tnf);
		}
	}
}

void GrSimulation::adjustFauxDegradation(double dt, bool tnfrDynamics, bool il10rDynamics)
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
		case DIFF_REC_EQ:
			_pDiffusion = new GrDiffusionFTCS();
			break;
		case DIFF_SOR_CORRECT:
			_pDiffusion = new GrDiffusionBTCS();
			break;
		case DIFF_SOR_WRONG:
			_pDiffusion = new GrDiffusionWrongBTCS();
			break;
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
