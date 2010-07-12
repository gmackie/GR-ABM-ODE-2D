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
#include "areatest.h"
#include "mtbtest.h"
#include "recruitmentprob.h"

GrSimulation::GrSimulation()
	: _time(0)
	, _grid()
	, _macList()
	, _tgamList()
	, _tcytList()
	, _tregList()
	, _stats()
	, _areaThreshold(DBL_MAX)
	, _pDiffusion(new GrDiffusionWrongBTCS())
	, _pTTest()
	, _pRecruitment(NULL)
	, _tnfrDynamics(false)
	
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

	// serialize time
	out << _time << std::endl;

	// serialize area threshold
	out << _areaThreshold << std::endl;

	// serialize diffusion method
	int intVal = (int) _pDiffusion->getMethod();
	out << intVal << std::endl;

	// serialize random number generator
	g_Rand.serialize(out);

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
}

void GrSimulation::deserialize(std::istream& in)
{
	assert(in.good());

	int intVal;

	// deserialize time
	in >> _time;

	// deserialize area threshold
	in >> _areaThreshold;

	// deserialize diffusion method
	in >> intVal;
	setDiffusionMethod((DiffusionMethod) intVal);

	// deserialize random number generator
	g_Rand.deserialize(in);

	// deserialize grid
	_grid.deserialize(in);

	// deserialize macs
	_macList.clear();
	in >> intVal;
	for (int i = 0; i < intVal; i++)
	{
		// create a dummy mac
		_macList.push_back(Mac(0, 0, 0, MAC_DEAD, 0, false, false));
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
		_tgamList.push_back(Tgam(0, 0, 0, TGAM_DEAD));
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
		_tcytList.push_back(Tcyt(0, 0, 0, TCYT_DEAD));
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
		_tregList.push_back(Treg(0, 0, 0, TREG_DEAD));
		Treg* pTreg = &_tregList.back();

		// update attributes of dummy mac and add to grid
		pTreg->deserialize(in);
		assert_res(_grid(pTreg->getRow(), pTreg->getCol()).addAgent(pTreg));
	}

	// deserialize statistics
	_stats.deserialize(in);
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
	while (count > 0)
	{
		int row = g_Rand.getInt(NROWS);
		int col = g_Rand.getInt(NCOLS);
		if (!_grid(row, col).hasMac())
		{
			createMac(row, col, g_Rand.getInt(-1, -maxMacAge), MAC_RESTING, false, false);
			count--;
		}
	}
}

void GrSimulation::solve()
{
	// here we perform a timestep, which is 10 minutes
	_time++;

	// reset statistics
	_stats.reset();

	// calculate diffusion every 10 minutes
	// dt = 6s, solve for 10 minutes = 100 * 6 seconds
	double dt = 6;
	for (int t = 0; t < 100; t++) 
	{
		secreteFromMacrophages();
		secreteFromCaseations();
		_pDiffusion->diffuse(_grid);
		if (_tnfrDynamics)
		{
			updateReceptorDynamics(dt);
		}
		else
		{
			adjustTNFDegradation(dt);
		}
	}
	
	// move macrophages
	moveMacrophages();

	// move T cells every 10 minutes
	moveTcells();
	
	// recruit agents from vascular sources every 10 minutes
	_pRecruitment->recruit(*this);

	// reset statistics
	_stats.resetAgentStats();

	// compute next state every 10 minutes
	computeNextStates();

	// update extracellular Mtb
	growExtMtb();

	// update states and remove dead agents from lists and grid
	updateStates();

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
		_stats.updateMacStatistics(mac.getState());

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
		_stats.updateTgamStatistics(tgam.getState());

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
		_stats.updateTcytStatistics(tcyt.getState());

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
		_stats.updateTregStatistics(treg.getState());

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
		it->computeNextState(_time, _grid, _stats, _tnfrDynamics);
	}
	for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
	{
		it->computeNextState(_time, _grid, _stats, _tnfrDynamics);
	}
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
		it->computeNextState(_time, _grid, _stats, _tnfrDynamics);
	}
	for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
	{
		it->computeNextState(_time, _grid, _stats, _tnfrDynamics);
	}
}

void GrSimulation::secreteFromMacrophages()
{
	for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
	{
		it->secrete(_grid);
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
			}
		}
	}
}

void GrSimulation::updateReceptorDynamics(double dt)
{
	for (MacList::iterator it = _macList.begin(); it != _macList.end(); it++)
	{
		it->solveODEs(_grid, dt);
	}
	for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
	{
		it->solveODEs(_grid, dt);
	}
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
		it->solveODEs(_grid, dt);
	}
	for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
	{
		it->solveODEs(_grid, dt);
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
			it->move(_grid);
		}
		// move active macrophages
		else if (_time % timeActive == 0 && it->getState() == MAC_ACTIVE)
		{
			it->move(_grid);
		}
		// move infected macrophages
		else if (_time % timeInfected == 0 && it->getState() == MAC_INFECTED)
		{
			it->move(_grid);
		}
	}
}

void GrSimulation::moveTcells()
{
	for (TgamList::iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
	{
		it->move(_grid);
	}
	for (TcytList::iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
		it->move(_grid);
	}
	for (TregList::iterator it = _tregList.begin(); it != _tregList.end(); it++)
	{
		it->move(_grid);
	}
}

void GrSimulation::growExtMtb()
{
	const double growthRate = _PARAM(PARAM_EXTMTB_GROWTH_RATE) - 1;
	const double upperBound = _PARAM(PARAM_EXTMTB_UPPER_BOUND) * 1.1;

	for (int row = 0; row < NROWS; row++)
	{
		for (int col = 0; col < NCOLS; col++)
		{
			GridCell& cell = _grid(row, col);
			
			double extMtb = cell.getExtMtb();
			double dExtMtb = growthRate * extMtb * (1.0 - extMtb / upperBound);
			
			cell.incExtMtb(dExtMtb);

			_stats.incTotExtMtb(cell.getExtMtb());
			_stats.incTotMacAttractant(cell.getMacAttractant());
			_stats.incTotTNF(cell.getTNF());
			_stats.incTotCCL2(cell.getCCL2());
			_stats.incTotCCL5(cell.getCCL5());
			_stats.incTotCXCL9(cell.getCXCL9());

			if (cell.getTNF() >= _areaThreshold)
				_stats.incArea();
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
			
			// simulate the effect of TNF internalization by cells in the form of degradation
			double dtnf;
			double tnf = cell.getTNF();
			if (cell.hasMac())
			{
				dtnf = -_PARAM(PARAM_GR_K_INT1) * (tnf / (tnf + _PARAM(PARAM_GR_KD1) * 48.16e11)) * 1500 * dt * 0.4;
				tnf += dtnf;
			}
			dtnf = -_PARAM(PARAM_GR_K_INT1) * (tnf / (tnf + _PARAM(PARAM_GR_KD1) * 48.16e11)) * 800 * (cell.hasTcell()) * dt * 0.4;
			tnf += dtnf;
			
			cell.setTNF(tnf);
		}
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
