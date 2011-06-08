/*
 * tgamma.cpp
 *
 *  Created on: 03-nov-2009
 *      Author: M. El-Kebir
 */

#include "tgamma.h"
#include "grgrid.h"

// Needed for deserializing the model state.
// Avoids the calls to the random number generator in the normal constructor, allowing the random number generator
// to remain in synch after deserialization.
Tgam::Tgam()
	: Tcell()
	, _state(TGAM_DEAD)
	, _nextState(TGAM_DEAD)
	, _deactivationTime(-1)
	, _mTNF(-1.0)
	, _surfTNFR1(-1.0)
	, _surfTNFR2(-1.0)
	, _surfBoundTNFR1(-1.0)
	, _surfBoundTNFR2(-1.0)
	, _intBoundTNFR1(-1.0)
	, _intBoundTNFR2(-1.0)
	, _vTNFR1(-1.0)
	, _vTNFR2(-1.0)
	, _kSynth(-1.0)
	, _kTACE(-1.0)
{
}

Tgam::Tgam(int birthtime, int row, int col, TgamState state)
	: Tcell(birthtime, row, col)
	, _state(state)
	, _nextState(state)
	, _deactivationTime(-1)
	, _mTNF(0.0)
	, _surfTNFR1(g_Rand.getReal(_PARAM(PARAM_GR_MIN_TNFR1_TCELL),_PARAM(PARAM_GR_MAX_TNFR1_TCELL)))
	, _surfTNFR2(g_Rand.getReal(_PARAM(PARAM_GR_MIN_TNFR2_TCELL),_PARAM(PARAM_GR_MAX_TNFR2_TCELL)))
	, _surfBoundTNFR1(0.0)
	, _surfBoundTNFR2(0.0)
	, _intBoundTNFR1(0.0)
	, _intBoundTNFR2(0.0)
	, _vTNFR1(_surfTNFR1 * _PARAM(PARAM_GR_K_T1))
	, _vTNFR2(_surfTNFR2 * _PARAM(PARAM_GR_K_T2))
	, _kSynth(_PARAM(PARAM_GR_K_SYNTH_TCELL))
	, _kTACE(_PARAM(PARAM_GR_K_TACE_TCELL))
{
}

Tgam::~Tgam()
{
}

void Tgam::move(GrGrid& grid)
{
	Tcell::moveTcell(grid, true, true, true);
}

void Tgam::secrete(GrGrid& grid, bool tnfrDynamics, bool nfkbDynamics, bool tnfKnockout)
{
	if (_deactivationTime != -1)
	{
		_kSynth = 0;
		return;
	}
	
	GridCell& cell = grid(_row, _col);
	_kSynth = _PARAM(PARAM_GR_K_SYNTH_TCELL);
	if (!tnfrDynamics && !tnfKnockout)
		cell.incTNF(_PARAM(PARAM_TGAM_SEC_RATE_TNF));
}

void Tgam::computeNextState(const int time, GrGrid& grid, GrStat& stats, bool tnfrDynamics, bool nfkbDynamics)
{
	GridCell& cell = grid(_row, _col);
	double tnfBoundFraction = cell.getTNF() / (cell.getTNF() + _PARAM(PARAM_GR_KD1) * 48.16e11);

	// check if it is time to die
	if (timeToDie(time))
	{
		_nextState = TGAM_DEAD;
	}
	else if (tnfrDynamics && _intBoundTNFR1 > _PARAM(PARAM_GR_THRESHOLD_APOPTOSIS_TNF_MOLECULAR) &&	
			 g_Rand.getReal() < 1 - pow(2.7183, -_PARAM(PARAM_GR_K_APOPTOSIS_MOLECULAR) * (_intBoundTNFR1 - _PARAM(PARAM_GR_THRESHOLD_APOPTOSIS_TNF_MOLECULAR))))
	{
		// TNF induced apoptosis
		_nextState = TGAM_DEAD;
	}
	else if (!tnfrDynamics && tnfBoundFraction > _PARAM(PARAM_GR_THRESHOLD_APOPTOSIS_TNF) &&
			 g_Rand.getReal() < 1 - pow(2.7183, -_PARAM(PARAM_GR_K_APOPTOSIS) * (tnfBoundFraction - _PARAM(PARAM_GR_THRESHOLD_APOPTOSIS_TNF))))
	{
		// TNF induced apoptosis
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
	
	if (cell.hasMac())
	{
		// get the macrophage
		Mac* pMac = dynamic_cast<Mac*>(cell.getAgent(0));
		if (!pMac) pMac = dynamic_cast<Mac*>(cell.getAgent(1));

		// If the mac died on this time step ignore it.
		if (pMac->getNextState() == MAC_DEAD)
		{
			return;
		}

		_nextState = TGAM_ACTIVE;

		// Fas/FasL induced apoptosis with probability
		if (pMac &&
			(pMac->getState() == MAC_INFECTED || pMac->getState() == MAC_CINFECTED) &&
			g_Rand.getReal() < _PARAM(PARAM_TGAM_PROB_APOPTOSIS_FAS_FASL))
		{
			stats.incApoptosisFasFasL();
			pMac->apoptosis(grid);
			pMac->kill();

			cell.incNrKillings();
		}
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

void Tgam::solveODEs(GrGrid& grid, double dt)
{
	GridCell& cell = grid(_row, _col);
	
	double koff1 = _PARAM(PARAM_GR_K_ON1) * _PARAM(PARAM_GR_KD1);
	double koff2 = _PARAM(PARAM_GR_K_ON2) * _PARAM(PARAM_GR_KD2);
	double density = 1.25e11; // used for conversion of conc. unit (M -> #/cell) based on cell and microcompartment volumes 
	double Nav = 6.02e23; // Avogadro Number
	double vol = 8.0e-12; // volume of a cell in liter
	
	double tnf = cell.getTNF() / (Nav * vol);
	double shedtnfr2 = cell.getShedTNFR2() / (Nav * vol);
	
	double dmTNF;
	double dsurfTNFR1;
	double dsurfTNFR2;
	double dsurfBoundTNFR1;
	double dsurfBoundTNFR2;
	double dintBoundTNFR1;
	double dintBoundTNFR2;
	double dsTNF;
	double dshedTNFR2;
	
	dmTNF = (_kSynth - _kTACE * _mTNF) * dt;
	dsurfTNFR1 = (_vTNFR1 - _PARAM(PARAM_GR_K_ON1) * tnf * _surfTNFR1 + koff1 * _surfBoundTNFR1 - _PARAM(PARAM_GR_K_T1) * _surfTNFR1 + _PARAM(PARAM_GR_K_REC1) * _intBoundTNFR1) * dt;
	dsurfTNFR2 = (_vTNFR2 - _PARAM(PARAM_GR_K_ON2) * tnf * _surfTNFR2 + koff2 * _surfBoundTNFR2 - _PARAM(PARAM_GR_K_T2) * _surfTNFR2 + _PARAM(PARAM_GR_K_REC2) * _intBoundTNFR2) * dt;
	dsurfBoundTNFR1 = (_PARAM(PARAM_GR_K_ON1) * tnf * _surfTNFR1 - koff1 * _surfBoundTNFR1 - _PARAM(PARAM_GR_K_INT1) * _surfBoundTNFR1) * dt;
	dsurfBoundTNFR2 = (_PARAM(PARAM_GR_K_ON2) * tnf * _surfTNFR2 - koff2 * _surfBoundTNFR2 - _PARAM(PARAM_GR_K_INT2) * _surfBoundTNFR2 - _PARAM(PARAM_GR_K_SHED) * _surfBoundTNFR2) * dt;
	dintBoundTNFR1 = (_PARAM(PARAM_GR_K_INT1) * _surfBoundTNFR1 - _PARAM(PARAM_GR_K_DEG1) * _intBoundTNFR1 - _PARAM(PARAM_GR_K_REC1) * _intBoundTNFR1) * dt;
	dintBoundTNFR2 = (_PARAM(PARAM_GR_K_INT2) * _surfBoundTNFR2 - _PARAM(PARAM_GR_K_DEG2) * _intBoundTNFR2 - _PARAM(PARAM_GR_K_REC2) * _intBoundTNFR2) * dt;
	dsTNF = ((density/Nav) * (_kTACE * _mTNF - _PARAM(PARAM_GR_K_ON1) * tnf * _surfTNFR1 + koff1 * _surfBoundTNFR1 - _PARAM(PARAM_GR_K_ON2) * tnf * _surfTNFR2 + koff2 * _surfBoundTNFR2)) * dt; 
	dshedTNFR2 = ((density/Nav) * _PARAM(PARAM_GR_K_SHED) * _surfBoundTNFR2) * dt;
	
	_mTNF += dmTNF;
	_surfTNFR1 += dsurfTNFR1;
	_surfTNFR2 += dsurfTNFR2;
	_surfBoundTNFR1 += dsurfBoundTNFR1;
	_surfBoundTNFR2 += dsurfBoundTNFR2;
	_intBoundTNFR1 += dintBoundTNFR1;
	_intBoundTNFR2 += dintBoundTNFR2;
	tnf += dsTNF;
	shedtnfr2 += dshedTNFR2;
	
	cell.setTNF(Nav * vol * tnf);
	cell.setShedTNFR2(Nav * vol * shedtnfr2);
	if (_mTNF < 0 || _surfTNFR1 < 0 || _surfBoundTNFR1 < 0 || _surfTNFR2 < 0 || _surfBoundTNFR2 < 0)
		std::cout << "Error: Negative Value of Species in TNF/TNFR dynamics" << std::endl;
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
	out << _mTNF << std::endl;
	out << _surfTNFR1 << std::endl;
	out << _surfTNFR2 << std::endl;
	out << _surfBoundTNFR1 << std::endl;
	out << _surfBoundTNFR2 << std::endl;
	out << _intBoundTNFR1 << std::endl;
	out << _intBoundTNFR2 << std::endl;
	out << _vTNFR1 << std::endl;
	out << _vTNFR2 << std::endl;
	out << _kSynth << std::endl;
	out << _kTACE << std::endl;
	
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
	in >> _mTNF;
	in >> _surfTNFR1;
	in >> _surfTNFR2;
	in >> _surfBoundTNFR1;
	in >> _surfBoundTNFR2;
	in >> _intBoundTNFR1;
	in >> _intBoundTNFR2;
	in >> _vTNFR1;
	in >> _vTNFR2;
	in >> _kSynth;
	in >> _kTACE;
}
