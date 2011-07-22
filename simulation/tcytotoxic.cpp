/*
 * tcytotoxic.cpp
 *
 *  Created on: 05-nov-2009
 *      Author: M. El-Kebir
 */

#include "tcytotoxic.h"
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

Tcyt::Tcyt(int birthtime, int row, int col, TcytState state)
	: Tcell(birthtime, row, col)
	, _state(state)
	, _nextState(state)
	, _deactivationTime(-1)
	, _mTNF(0.0)
	, _surfTNFR1(g_Rand.getReal(_PARAM(PARAM_GR_MEAN_TNFR1_TCELL)-_PARAM(PARAM_GR_STD_TNFR1_TCELL),_PARAM(PARAM_GR_MEAN_TNFR1_TCELL)+_PARAM(PARAM_GR_STD_TNFR1_TCELL)))
	, _surfTNFR2(g_Rand.getReal(_PARAM(PARAM_GR_MEAN_TNFR2_TCELL)-_PARAM(PARAM_GR_STD_TNFR2_TCELL),_PARAM(PARAM_GR_MEAN_TNFR2_TCELL)+_PARAM(PARAM_GR_STD_TNFR2_TCELL)))
//	, _surfTNFR1(g_Rand.getLogNormal(_PARAM(PARAM_GR_MEAN_TNFR1_TCELL),_PARAM(PARAM_GR_STD_TNFR1_TCELL)))
//	, _surfTNFR2(g_Rand.getLogNormal(_PARAM(PARAM_GR_MEAN_TNFR2_TCELL),_PARAM(PARAM_GR_STD_TNFR2_TCELL)))
	, _surfBoundTNFR1(0.0)
	, _surfBoundTNFR2(0.0)
	, _intBoundTNFR1(0.0)
	, _intBoundTNFR2(0.0)
	, _vTNFR1(_surfTNFR1 * _PARAM(PARAM_GR_K_T1))
	, _vTNFR2(_surfTNFR2 * _PARAM(PARAM_GR_K_T2))
	, _kSynth(_PARAM(PARAM_GR_K_SYNTH_TCELL)/10)
	, _kTACE(_PARAM(PARAM_GR_K_TACE_TCELL))
{
}

Tcyt::~Tcyt()
{
}

void Tcyt::move(GrGrid& grid)
{
	Tcell::moveTcell(grid, false, true, true);
}

void Tcyt::secrete(GrGrid& grid, bool tnfrDynamics, bool nfkbDynamics, bool tnfDepletion)
{
	if (_deactivationTime != -1)
	{
		_kSynth = 0;
		return;
	}
	
	GridCell& cell = grid(_row, _col);
	_kSynth = _PARAM(PARAM_GR_K_SYNTH_TCELL)/10;
	if (!tnfrDynamics && !tnfDepletion)
		cell.incTNF(_PARAM(PARAM_TCYT_SEC_RATE_TNF));
}

void Tcyt::computeNextState(const int time, GrGrid& grid, GrStat& stats, bool tnfrDynamics, bool nfkbDynamics)
{
	GridCell& cell = grid(_row, _col);
	double tnfBoundFraction = cell.getTNF() / (cell.getTNF() + _PARAM(PARAM_GR_KD1) * 48.16e11);

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
					// kill, half the intracellular bacteria disperse to the Moore neighborhood
					pMac->disperseMtb(grid, 0.5);

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

void Tcyt::solveODEs(GrGrid& grid, double dt)
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

	Serialization::writeHeader(out, Tcyt::_ClassName);

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

	if (!Serialization::readFooter(in, Tcyt::_ClassName))
	{
		exit(1);
	}
}
