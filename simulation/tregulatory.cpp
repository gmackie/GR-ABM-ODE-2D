/*
 * tregulatory.cpp
 *
 *  Created on: 05-nov-2009
 *      Author: M. El-Kebir
 */

#include "tregulatory.h"
#include "grgrid.h"
#include "serialization.h"

const std::string Treg::_ClassName = "Treg";


// Needed for deserializing the model state.
// Avoids the calls to the random number generator in the normal constructor, allowing the random number generator
// to remain in synch after deserialization.
Treg::Treg()
	: Tcell()
	, _state(TREG_DEAD)
	, _nextState(TREG_DEAD)
{
}

Treg::Treg(int birthtime, int row, int col, TregState state)

	: Tcell(birthtime, row, col, 0.0 /* kSynth */)
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

void Treg::secrete(GrGrid& grid, bool, bool, bool, bool il10rDynamics, bool il10Depletion)
{
    
    _kISynth = _PARAM(PARAM_GR_I_K_SYNTH_TCELL);
    
    if (!il10rDynamics && !il10Depletion) {
        grid.incil10(_pos, (_PARAM(PARAM_TREG_SEC_RATE_IL10)));
    }
    
}

void Treg::deactivate(const int)
{
}

void Treg::computeNextState(const int time, GrGrid& grid, GrStat& stats, bool tnfrDynamics, bool, bool, bool)
{
	double tnfBoundFraction = grid.TNF(_pos) / (grid.TNF(_pos) + _PARAM(PARAM_GR_KD1) * 48.16e11);

	// check if it is time to die
	if (timeToDie(time))
	{
		_nextState = TREG_DEAD;
	}
	else if (tnfrDynamics && _intBoundTNFR1 > _PARAM(PARAM_GR_THRESHOLD_APOPTOSIS_TNF_MOLECULAR) &&	
			 g_Rand.getReal() < 1 - pow(2.7183, -_PARAM(PARAM_GR_K_APOPTOSIS_MOLECULAR) * (_intBoundTNFR1 - _PARAM(PARAM_GR_THRESHOLD_APOPTOSIS_TNF_MOLECULAR))))
	{
		// TNF induced apoptosis
		stats.incTcellApoptosisTNF();
		_nextState = TREG_DEAD;
	}
	else if (!tnfrDynamics && tnfBoundFraction > _PARAM(PARAM_GR_THRESHOLD_APOPTOSIS_TNF) &&
			 g_Rand.getReal() < 1 - pow(2.7183, -_PARAM(PARAM_GR_K_APOPTOSIS) * (tnfBoundFraction - _PARAM(PARAM_GR_THRESHOLD_APOPTOSIS_TNF))))
	{
		// TNF induced apoptosis
		stats.incTcellApoptosisTNF();
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
    default: throw std::runtime_error("Unknown Tcyt state");
		}
	}
}

void Treg::handleResting(const int time, GrGrid& grid, GrStat&)
{
	for (int i = -1; i <= 1; i++)
	{
		for (int j = -1; j <= 1; j++)
		{
      Pos p(grid.mod_row(_pos.x+i), grid.mod_col(_pos.y+j));
			Agent* pAgent0 = grid.agent(p, 0);
			Agent* pAgent1 = grid.agent(p, 1);

			if (pAgent0 && g_Rand.getReal() <= _PARAM(PARAM_TREG_PROB_DOWN_REGULATE))
			{
				pAgent0->deactivate(time);
			}
			else if (pAgent1 && g_Rand.getReal() <= _PARAM(PARAM_TREG_PROB_DOWN_REGULATE))
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
  default: throw std::runtime_error("Unknown Tcyt state");
	}
}

void Treg::serialize(std::ostream& out) const
{
	assert(out.good());

	Serialization::writeHeader(out, Treg::_ClassName);

	Tcell::serialize(out);

	int intVal = (int) _state;
	out << intVal << std::endl;

	intVal = (int) _nextState;
	out << intVal << std::endl;

	Serialization::writeFooter(out, Treg::_ClassName);
}

void Treg::deserialize(std::istream& in)
{
	assert(in.good());

	if (!Serialization::readHeader(in, Treg::_ClassName))
	{
		exit(1);
	}

	int intVal;

	Tcell::deserialize(in);

	in >> intVal;
	_state = (TregState) intVal;

	in >> intVal;
	_nextState = (TregState) intVal;

	if (!Serialization::readFooter(in, Treg::_ClassName))
	{
		exit(1);
	}
}
