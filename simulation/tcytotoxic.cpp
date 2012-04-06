/*
 * tcytotoxic.cpp
 *
 *  Created on: 05-nov-2009
 *      Author: M. El-Kebir
 */

#include "tcytotoxic.h"
#include "macrophage.h"
#include "grgrid.h"
#include "serialization.h"

const std::string Tcyt::_ClassName = "Tcyt";

// Needed for deserializing the model state.
// Avoids the calls to the random number generator in the normal constructor, allowing the random number generator
// to remain in synch after deserialization.
Tcyt::Tcyt()
	: Tcell()
	, _state(TCYT_DEAD)
	, _nextState(TCYT_DEAD)
	, _deactivationTime(-1)
{
}

Tcyt::Tcyt(int birthtime, int row, int col, TcytState state)

	: Tcell(birthtime, row, col, _PARAM(PARAM_GR_K_SYNTH_TCELL)/10)
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

void Tcyt::secrete(GrGrid& grid, bool tnfrDynamics, bool, bool tnfDepletion, bool il10rDynamics, bool il10Depletion, int mdt)
{
	if (_deactivationTime != -1)
	{
		_kSynth = 0;
        _kmRNA = 0;
        _kISynth = 0;
		return;
	}
	
	_kSynth = _PARAM(PARAM_GR_K_SYNTH_TCELL)/10;
    _kmRNA = _PARAM(PARAM_GR_K_RNA_TCELL)/10;
    _kISynth = 0.0;
    
	if (!tnfrDynamics && !tnfDepletion)
    {    
       const double il10 = log(((grid.il10(_pos) * MW_IL10 * 1e6)/(NAV * VOL))); // converting il10 concentration to log(ng/mL) for use in dose dependence
       const double tnfMOD = (1.0/(1.0 + exp((il10 + _PARAM(PARAM_GR_LINK_LOG_ALPHA))/_PARAM(PARAM_GR_LINK_LOG_BETA)))); // calculate the fraction of inhibition
        
		grid.incTNF(_pos, (tnfMOD * _PARAM(PARAM_TCYT_SEC_RATE_TNF) * mdt));
    }
    if (!il10rDynamics && !il10Depletion) {
        grid.setil10(_pos, (_PARAM(PARAM_TCYT_SEC_RATE_IL10) * mdt));
    }
    
}

void Tcyt::computeNextState(const int time, GrGrid& grid, GrStat& stats, bool tnfrDynamics, bool, bool, bool)
{
	double tnfBoundFraction = grid.TNF(_pos) / (grid.TNF(_pos) + _PARAM(PARAM_GR_KD1) * 48.16e11);

	// check if it is time to die
	if (timeToDie(time))
	{
		_nextState = TCYT_DEAD;
	}
	else if (tnfrDynamics && _intBoundTNFR1 > _PARAM(PARAM_GR_THRESHOLD_APOPTOSIS_TNF_MOLECULAR) &&	
			 g_Rand.getReal() < 1 - pow(2.7183, -_PARAM(PARAM_GR_K_APOPTOSIS_MOLECULAR) * (_intBoundTNFR1 - _PARAM(PARAM_GR_THRESHOLD_APOPTOSIS_TNF_MOLECULAR))))
	{
		// TNF induced apoptosis
		stats.incTcellApoptosisTNF();
		_nextState = TCYT_DEAD;
	}
	else if (!tnfrDynamics && tnfBoundFraction > _PARAM(PARAM_GR_THRESHOLD_APOPTOSIS_TNF) &&
			 g_Rand.getReal() < 1 - pow(2.7183, -_PARAM(PARAM_GR_K_APOPTOSIS) * (tnfBoundFraction - _PARAM(PARAM_GR_THRESHOLD_APOPTOSIS_TNF))))
	{
		// TNF induced apoptosis
		stats.incTcellApoptosisTNF();
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
    default:
      throw std::runtime_error("Unknown Tcyt state");
		}
	}
}

void Tcyt::handleActive(const int, GrGrid& grid, GrStat&)
{

	if (grid.hasAgentType(MAC, _pos))
	{
		Mac* pMac = dynamic_cast<Mac*>(grid.agent(_pos, 0));
		if (!pMac) pMac = dynamic_cast<Mac*>(grid.agent(_pos, 1));

		assert(pMac);

		// If the mac died on this time step ignore it.
		if (pMac->getNextState() == MAC_DEAD)
		{
			return;
		}

		if (g_Rand.getReal() < _PARAM(PARAM_TCYT_PROB_KILL_MAC))
		{
			if (pMac->getState() == MAC_INFECTED)
			{
				pMac->setIntMtb(0);
				pMac->kill();

				// contribute to caseation
				if (!grid.incKillings(_pos))
					_nextState = TCYT_ACTIVE;
			}
			else if (pMac->getState() == MAC_CINFECTED)
			{
				double r = g_Rand.getReal();
				if (r < _PARAM(PARAM_TCYT_PROB_KILL_MAC_CLEANLY))
				{
					pMac->setIntMtb(0);
					pMac->kill();
					if (!grid.incKillings(_pos))
						_nextState = TCYT_ACTIVE;
				}
				else
				{
					// kill, half the intracellular bacteria disperse to the Moore neighborhood
					pMac->disperseMtb(grid, 0.5);

					pMac->setIntMtb(0);
					pMac->kill();
					if (!grid.incKillings(_pos))
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
  default: throw std::runtime_error("Unknown Tcyt state");
	}
}

void Tcyt::serialize(std::ostream& out) const
{
	assert(out.good());

	Serialization::writeHeader(out, Tcyt::_ClassName);

	Tcell::serialize(out);

	int intVal = (int) _state;
	out << intVal << std::endl;

	intVal = (int) _nextState;
	out << intVal << std::endl;

	out << _deactivationTime << std::endl;

	Serialization::writeFooter(out, Tcyt::_ClassName);
}

void Tcyt::deserialize(std::istream& in)
{
	assert(in.good());

	if (!Serialization::readHeader(in, Tcyt::_ClassName))
	{
		exit(1);
	}

	int intVal;

	Tcell::deserialize(in);

	in >> intVal;
	_state = (TcytState) intVal;

	in >> intVal;
	_nextState = (TcytState) intVal;

	in >> _deactivationTime;

	if (!Serialization::readFooter(in, Tcyt::_ClassName))
	{
		exit(1);
	}
}
