/*
 * tgamma.cpp
 *
 *  Created on: 03-nov-2009
 *      Author: M. El-Kebir
 */

#include "tgamma.h"
#include "grgrid.h"

Tgam::Tgam(int birthtime, int row, int col, TgamState state)
	: Tcell(birthtime, row, col)
	, _state(state)
	, _nextState(state)
	, _deactivationTime(-1)
{
}

Tgam::~Tgam()
{
}

void Tgam::move(GrGrid& grid)
{
	Tcell::moveTcell(grid, true, true, true);
}

void Tgam::secrete(GrGrid&)
{
	// a macrophage calls this function
	//grid(_row, _col).incTNF(Params::_INCR_TNF_TGAM);
	_state = TGAM_ACTIVE;
}

void Tgam::computeNextState(const int time, GrGrid& grid, GrStat& stats)
{
	GridCell& cell = grid(_row, _col);

	// check if it is time to die
	if (timeToDie(time))
	{
		_nextState = TGAM_DEAD;
	}
	else if (cell.getTNF() > _PARAM(PARAM_GR_THRESHOLD_APOPTOSIS_TNF) &&	
			 g_Rand.getReal() < 1 - pow(2.7183, -_PARAM(PARAM_GR_K_APOPTOSIS) * (cell.getTNF() - _PARAM(PARAM_GR_THRESHOLD_APOPTOSIS_TNF))))
	{
		// TNF induced apoptosis (with probability PARAM_GR_PROB_APOPTOSIS_TNF)
		_nextState = TGAM_DEAD;
	}
	else
	{
		switch (_state)
		{
		case TGAM_DEAD:
			// if dead, stay dead
			_nextState = TGAM_DEAD;
			break;
		case TGAM_ACTIVE:
			handleActive(time, grid, stats);
			break;
		case TGAM_DOWN_REGULATED:
			handleDownRegulated(time, grid, stats);
			break;
		}
	}
}

void Tgam::handleActive(const int, GrGrid& grid, GrStat& stats)
{
	GridCell& cell = grid(_row, _col);
	
	// get the macrophage
	Mac* pMac = dynamic_cast<Mac*>(cell.getAgent(0));
	if (!pMac) pMac = dynamic_cast<Mac*>(cell.getAgent(1));

	_nextState = TGAM_ACTIVE;

	// Fas/FasL induced apoptosis with probability 
	if (pMac &&
		(pMac->getState() == MAC_INFECTED || pMac->getState() == MAC_CINFECTED) &&
		g_Rand.getReal() < _PARAM(PARAM_TGAM_PROB_APOPTOSIS_FAS_FASL))
	{
		stats.incApoptosisFasFasL();
		pMac->apoptosis(grid);
		pMac->kill();

		/* Does Fas/FasL induced apoptosis contribute to caseation? */
		cell.incNrKillings();
	}
}

void Tgam::handleDownRegulated(const int time, GrGrid&, GrStat&)
{
	if (time - _deactivationTime >= _PARAM(PARAM_TGAM_TIMESPAN_REGULATED))
	{
		_nextState = TGAM_ACTIVE;
	}
	else
	{
		_nextState = TGAM_DOWN_REGULATED;
	}
}

void Tgam::deactivate(const int time)
{
	_nextState = _state = TGAM_DOWN_REGULATED;
	_deactivationTime = time;
}

void Tgam::updateState()
{
	_state = _nextState;
}

void Tgam::kill()
{
	_nextState = _state = TGAM_DEAD;
}

void Tgam::print() const
{
	std::cout << "Tgam - " << _birthTime << " - ";

	switch (_state)
	{
	case TGAM_DEAD:
		std::cout << "dead" << std::endl;
		break;
	case TGAM_DOWN_REGULATED:
		std::cout << "down-regulated" << std::endl;
		break;
	case TGAM_ACTIVE:
		std::cout << "active" << std::endl;
		break;
	}
}

void Tgam::serialize(std::ostream& out) const
{
	assert(out.good());

	Tcell::serialize(out);

	int intVal = (int) _state;
	out << intVal << std::endl;

	intVal = (int) _nextState;
	out << intVal << std::endl;

	out << _deactivationTime << std::endl;
}

void Tgam::deserialize(std::istream& in)
{
	assert(in.good());

	int intVal;

	Tcell::deserialize(in);

	in >> intVal;
	_state = (TgamState) intVal;

	in >> intVal;
	_nextState = (TgamState) intVal;

	in >> _deactivationTime;
}
