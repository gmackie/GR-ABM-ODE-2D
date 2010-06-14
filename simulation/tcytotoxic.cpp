/*
 * tcytotoxic.cpp
 *
 *  Created on: 05-nov-2009
 *      Author: M. El-Kebir
 */

#include "tcytotoxic.h"
#include "grgrid.h"

Tcyt::Tcyt(int birthtime, int row, int col, TcytState state)
	: Tcell(birthtime, row, col)
	, _state(state)
	, _nextState(state)
	, _deactivationTime(-1)
{
}

Tcyt::~Tcyt()
{
}

void Tcyt::move(GrGrid& grid)
{
	Tcell::moveTcell(grid, false, true, true);
}

void Tcyt::secrete(GrGrid&)
{
}

void Tcyt::computeNextState(const int time, GrGrid& grid, GrStat& stats)
{
	GridCell& cell = grid(_row, _col);

	// check if it is time to die
	if (timeToDie(time))
	{
		_nextState = TCYT_DEAD;
	}
	else if (cell.getTNF() > _PARAM(PARAM_GR_THRESHOLD_APOPTOSIS_TNF) &&
		g_Rand.getReal() < _PARAM(PARAM_GR_PROB_APOPTOSIS_TNF))
	{
		// TNF induced apoptosis (with probability PARAM_GR_PROB_APOPTOSIS_TNF)
		_nextState = TCYT_DEAD;
	}
	else
	{
		switch (_state)
		{
		case TCYT_DEAD:
			// if dead, stay dead
			_nextState = TCYT_DEAD;
			break;
		case TCYT_ACTIVE:
			handleActive(time, grid, stats);
			break;
		case TCYT_DOWN_REGULATED:
			handleDownRegulated(time, grid, stats);
			break;
		}
	}
}

void Tcyt::handleActive(const int, GrGrid& grid, GrStat&)
{
	GridCell& cell = grid(_row, _col);

	if (cell.hasMac())
	{
		Mac* pMac = dynamic_cast<Mac*>(cell.getAgent(0));
		if (!pMac) pMac = dynamic_cast<Mac*>(cell.getAgent(1));

		assert(pMac);

		if (g_Rand.getReal() < _PARAM(PARAM_TCYT_PROB_KILL_MAC))
		{
			if (pMac->getState() == MAC_INFECTED)
			{
				pMac->setIntMtb(0);
				pMac->kill();

				// contribute to caseation
				if (!cell.incNrKillings())
					_nextState = TCYT_ACTIVE;
			}
			else if (pMac->getState() == MAC_CINFECTED)
			{
				double r = g_Rand.getReal();
				if (r < _PARAM(PARAM_TCYT_PROB_KILL_MAC_CLEANLY))
				{
					pMac->setIntMtb(0);
					pMac->kill();
					if (!cell.incNrKillings())
						_nextState = TCYT_ACTIVE;
				}
				else
				{
					// kill, intracellular bacteria disperse to the Moore neighborhood
					const double dExtMtb = pMac->getIntMtb() / 9.0;

					/* is this ok? none of the bacteria are killed. */
					for (int i = -1; i <= 1; i++)
						for (int j = -1; j <= 1; j++)
							grid(MOD_ROW(_row + i), MOD_COL(_col + j)).incExtMtb(dExtMtb);

					pMac->setIntMtb(0);
					pMac->kill();
					if (!cell.incNrKillings())
						_nextState = TCYT_ACTIVE;
				}
			}
		}
	}
}

void Tcyt::handleDownRegulated(const int time, GrGrid&, GrStat&)
{
	if (time - _deactivationTime >= _PARAM(PARAM_TCYT_TIMESPAN_REGULATED))
	{
		_nextState = TCYT_ACTIVE;
	}
	else
	{
		_nextState = TCYT_DOWN_REGULATED;
	}
}

void Tcyt::deactivate(const int time)
{
	_nextState = _state = TCYT_DOWN_REGULATED;
	_deactivationTime = time;
}

void Tcyt::updateState()
{
	_state = _nextState;
}

void Tcyt::kill()
{
	_nextState = _state = TCYT_DEAD;
}

void Tcyt::print() const
{
	std::cout << "Tcyt - " << _birthTime << " - ";

	switch (_state)
	{
	case TCYT_DEAD:
		std::cout << "dead" << std::endl;
		break;
	case TCYT_ACTIVE:
		std::cout << "active" << std::endl;
		break;
	case TCYT_DOWN_REGULATED:
		std::cout << "down-regulated" << std::endl;
		break;
	}
}

void Tcyt::serialize(std::ostream& out) const
{
	assert(out.good());

	Tcell::serialize(out);

	int intVal = (int) _state;
	out << intVal << std::endl;

	intVal = (int) _nextState;
	out << intVal << std::endl;

	out << _deactivationTime << std::endl;
}

void Tcyt::deserialize(std::istream& in)
{
	assert(in.good());

	int intVal;

	Tcell::deserialize(in);

	in >> intVal;
	_state = (TcytState) intVal;

	in >> intVal;
	_nextState = (TcytState) intVal;

	in >> _deactivationTime;
}
