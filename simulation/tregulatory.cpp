/*
 * tregulatory.cpp
 *
 *  Created on: 05-nov-2009
 *      Author: M. El-Kebir
 */

#include "tregulatory.h"
#include "grgrid.h"

Treg::Treg(int birthtime, int row, int col, TregState state)
	: Tcell(birthtime, row, col)
	, _state(state)
	, _nextState(state)
{
}

Treg::~Treg()
{
}

void Treg::move(GrGrid& grid)
{
	Tcell::moveTcell(grid, false, true, false);
}

void Treg::secrete(GrGrid&)
{
}

void Treg::deactivate(const int)
{
}

void Treg::computeNextState(const int time, GrGrid& grid, GrStat& stats)
{
	GridCell& cell = grid(_row, _col);

	// check if it is time to die
	if (timeToDie(time))
	{
		_nextState = TREG_DEAD;
	}
	else if (cell.getTNF() > _PARAM(PARAM_GR_THRESHOLD_APOPTOSIS_TNF) &&	
			 g_Rand.getReal() < 1 - pow(2.7183, -_PARAM(PARAM_GR_K_APOPTOSIS) * (cell.getTNF() - _PARAM(PARAM_GR_THRESHOLD_APOPTOSIS_TNF))))
	{
		// TNF induced apoptosis (with probability _PROB_TNF_APOPTOSIS)
		_nextState = TREG_DEAD;
	}
	else
	{
		switch (_state)
		{
		case TREG_DEAD:
			// if dead, stay dead
			_nextState = TREG_DEAD;
			break;
		case TREG_ACTIVE:
			handleResting(time, grid, stats);
			break;
		}
	}
}

void Treg::handleResting(const int time, GrGrid& grid, GrStat&)
{
	for (int i = -1; i <= 1; i++)
	{
		for (int j = -1; j <= 1; j++)
		{
			GridCell& cell = grid(MOD_ROW(_row + i), MOD_COL(_col + j));
			Agent* pAgent0 = cell.getAgent(0);
			Agent* pAgent1 = cell.getAgent(1);

			if (pAgent0)
			{
				pAgent0->deactivate(time);
			}
			else if (pAgent1)
			{
				pAgent1->deactivate(time);
			}
		}
	}
}

void Treg::updateState()
{
	_state = _nextState;
}

void Treg::kill()
{
	_nextState = _state = TREG_DEAD;
}

void Treg::print() const
{
	std::cout << "Treg - " << _birthTime << " - ";

	switch (_state)
	{
	case TREG_DEAD:
		std::cout << "dead" << std::endl;
		break;
	case TREG_ACTIVE:
		std::cout << "active" << std::endl;
		break;
	}
}

void Treg::serialize(std::ostream& out) const
{
	assert(out.good());

	Tcell::serialize(out);

	int intVal = (int) _state;
	out << intVal << std::endl;

	intVal = (int) _nextState;
	out << intVal << std::endl;
}

void Treg::deserialize(std::istream& in)
{
	assert(in.good());

	int intVal;

	Tcell::deserialize(in);

	in >> intVal;
	_state = (TregState) intVal;

	in >> intVal;
	_nextState = (TregState) intVal;
}
