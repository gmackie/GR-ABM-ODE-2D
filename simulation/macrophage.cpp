/*
 * macrophage.cpp
 *
 *  Created on: 05-nov-2009
 *      Author: M. El-Kebir
 */

#include "macrophage.h"
#include "params.h"
#include "grgrid.h"

Mac::Mac(int birthtime, int row, int col, MacState state, double intMtb, bool NFkB, bool stat1)
	: Agent(birthtime, birthtime + _PARAM(PARAM_MAC_AGE), row, col)
	, _state(state)
	, _nextState(state)
	, _intMtb(intMtb)
	, _NFkB(NFkB)
	, _stat1(stat1)
	, _activationTime(-1)
	, _deactivationTime(-1)
	, _mTNF(0.0)
	, _surfTNFR1(g_Rand.getReal(_PARAM(PARAM_GR_MIN_TNFR1_MAC),_PARAM(PARAM_GR_MAX_TNFR1_MAC)))
	, _surfTNFR2(g_Rand.getReal(_PARAM(PARAM_GR_MIN_TNFR2_MAC),_PARAM(PARAM_GR_MAX_TNFR2_MAC)))
	, _surfBoundTNFR1(0.0)
	, _surfBoundTNFR2(0.0)
	, _intBoundTNFR1(0.0)
	, _intBoundTNFR2(0.0)
	, _vTNFR1(_surfTNFR1 * _PARAM(PARAM_GR_K_T1))
	, _vTNFR2(_surfTNFR2 * _PARAM(PARAM_GR_K_T2))
	, _kSynth(0.0)
	, _kTACE(_PARAM(PARAM_GR_K_TACE_MAC))
{
}

Mac::~Mac()
{
}

void Mac::move(GrGrid& grid)
{
	int k = Agent::moveAgent(grid, true, true, false, true, 1.5);

	/**
	 * 0 1 2
	 * 3 4 5
	 * 6 7 8
	 */
	int dRow = k / 3 - 1;
	int dCol = k % 3 - 1;
	int newRow = MOD_ROW(_row + dRow);
	int newCol = MOD_COL(_col + dCol);

	GridCell& cell = grid(_row, _col);
	GridCell& newCell = grid(newRow, newCol);
	if (newCell.addAgent(this))
	{
		_row = newRow;
		_col = newCol;
		assert_res(cell.removeAgent(this));
	}
}

void Mac::secrete(GrGrid& grid)
{
	if (_deactivationTime != -1)
		return;
	
	GridCell& cell = grid(_row, _col);
	
	if (_NFkB)
	{
		// if NFkB is turned on, secrete TNF and chemokines
		if (_state != MAC_RESTING)
		{
			cell.incCCL2(_PARAM(PARAM_MAC_SEC_RATE_CCL2));
			cell.incCCL5(_PARAM(PARAM_MAC_SEC_RATE_CCL5));
			cell.incCXCL9(_PARAM(PARAM_MAC_SEC_RATE_CXCL9));
			_kSynth = _PARAM(PARAM_GR_K_SYNTH_MAC);
			cell.incTNF(_PARAM(PARAM_MAC_SEC_RATE_TNF));
		}
		else
		{
			cell.incCCL2(0.5 * _PARAM(PARAM_MAC_SEC_RATE_CCL2));
			cell.incCCL5(0.5 * _PARAM(PARAM_MAC_SEC_RATE_CCL5));
			cell.incCXCL9(0.5 * _PARAM(PARAM_MAC_SEC_RATE_CXCL9));
			_kSynth = 0.5 * _PARAM(PARAM_GR_K_SYNTH_MAC);
			cell.incTNF(0.5 * _PARAM(PARAM_MAC_SEC_RATE_TNF));
		}
		
		cell.incNrSecretions();
	}
	else if (_state == MAC_RESTING)
		_kSynth = 0.0;
	else if (_state == MAC_INFECTED)
	{
		cell.incCCL2(0.5 * _PARAM(PARAM_MAC_SEC_RATE_CCL2));
		cell.incCCL5(0.5 * _PARAM(PARAM_MAC_SEC_RATE_CCL5));
		cell.incCXCL9(0.5 * _PARAM(PARAM_MAC_SEC_RATE_CXCL9));
		_kSynth = 0.5 * _PARAM(PARAM_GR_K_SYNTH_MAC);
		cell.incTNF(0.5 * _PARAM(PARAM_MAC_SEC_RATE_TNF));
		
		cell.incNrSecretions();
	}
}


void Mac::updateState()
{
	_state = _nextState;
}

void Mac::kill()
{
	_nextState = _state = MAC_DEAD;
}

void Mac::apoptosis(GrGrid& grid)
{
	_nextState = MAC_DEAD;

	/* Question: Does an infected macrophage undergoing apoptosis contribute to caseation? */
	if (_state == MAC_INFECTED || _state == MAC_CINFECTED)
	{
		// disperse bacteria to the Moore neighborhood
		// (half of the intracellular bacteria dies, the other half is dispersed
		const double dExtMtb = _intMtb / 18.0;

		for (int i = -1; i <= 1; i++)
			for (int j = -1; j <= 1; j++)
				grid(MOD_ROW(_row + i), MOD_COL(_col + j)).incExtMtb(dExtMtb);

		_intMtb = 0;
	}
}

void Mac::computeNextState(const int time, GrGrid& grid, GrStat& stats)
{
	GridCell& cell = grid(_row, _col);
	double tnfBoundFraction = cell.getTNF() / (cell.getTNF() + _PARAM(PARAM_GR_KD1) * 48.16e11);
	
	// check if it is time to die
	if (timeToDie(time))
	{
		// release intracellular bacteria
		cell.incExtMtb(_intMtb);
		_intMtb = 0;

		/* Should we call this caseation or necrosis? */
		// in case active, death contributes to caseation
		if (_state == MAC_ACTIVE)
		{
			cell.incNrKillings();
		}

		_nextState = MAC_DEAD;
	}
	else if (tnfBoundFraction > _PARAM(PARAM_GR_THRESHOLD_APOPTOSIS_TNF) &&
			 g_Rand.getReal() < 1 - pow(2.7183, -_PARAM(PARAM_GR_K_APOPTOSIS) * (tnfBoundFraction - _PARAM(PARAM_GR_THRESHOLD_APOPTOSIS_TNF))))
	{
		// TNF induced apoptosis
		stats.incApoptosisTNF();
		apoptosis(grid);
	}
	else
	{
		if (_deactivationTime != -1 &&
			_deactivationTime + _PARAM(PARAM_MAC_TIMESPAN_REGULATED) <= time)
		{
			_deactivationTime = -1;
		}

		if (_deactivationTime == -1)
		{
			// update _NFkB
			bool tnfInducedNFkB = tnfBoundFraction > _PARAM(PARAM_MAC_THRESHOLD_NFKB_TNF) && 
			g_Rand.getReal() < 1 - pow(2.7183, -_PARAM(PARAM_MAC_K_NFKB) * (tnfBoundFraction - _PARAM(PARAM_MAC_THRESHOLD_NFKB_TNF)));
			
			_NFkB = _state == MAC_CINFECTED || _state == MAC_ACTIVE || tnfInducedNFkB ||
				getExtMtbInMoore(grid) > _PARAM(PARAM_MAC_THRESHOLD_NFKB_EXTMTB);
				//cell.getExtMtb() > _PARAM(PARAM_MAC_THRESHOLD_NFKB_EXTMTB);

			switch (_state)
			{
			case MAC_DEAD:
				// if dead, stay dead
				_nextState = MAC_DEAD;
				break;
			case MAC_RESTING:
				handleResting(time, grid, stats);
				break;
			case MAC_INFECTED:
				handleInfected(time, grid, stats);
				break;
			case MAC_CINFECTED:
				handleChronicallyInfected(time, grid, stats);
				break;
			case MAC_ACTIVE:
				handleActivated(time, grid, stats);
				break;
			}
		}
	}
}

void Mac::solveODEs(GrGrid& grid, double dt)
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

void Mac::handleResting(const int time, GrGrid& grid, GrStat&)
{
	GridCell& cell = grid(_row, _col);

	_stat1 |= g_Rand.getReal() < getCountTgam(TGAM_ACTIVE, grid) * _PARAM(PARAM_MAC_PROB_STAT1_TGAM);

	/* You can get infected from only 1 bacteria;
	 * there should always be a probability associated with getting infected.
	 *
	 * Suggestion:
	 * - always engulf as many bacteria as specified by _EXT_MTB_ENGULF_MAC_REST
	 *   (another suggestion: we can even fix this parameter to 1 and thus get rid of it)
	 * - with a certain probability kill those bacteria, otherwise become infected */
	if (cell.getExtMtb() <= _PARAM(PARAM_MAC_NR_UPTAKE_RI_EXTMTB))
	{
		// kill extracellular bacteria, only if there are not too many
		cell.setExtMtb(0);
		_nextState = MAC_RESTING;
	}
	else if (g_Rand.getReal() < _PARAM(PARAM_MAC_PROB_KILL_R_EXTMTB))
	{
		// there are too many extracellular bacteria, but we can kill some of those
		// with probability PARAM_MAC_PROB_KILL_R_EXTMTB
		cell.incExtMtb(-1 * std::min(_PARAM(PARAM_MAC_NR_UPTAKE_RI_EXTMTB), cell.getExtMtb()));
		_nextState = MAC_RESTING;
	}
	else
	{
		// too many bacteria and no killing => macrophage may become infected
		_intMtb = _PARAM(PARAM_MAC_NR_UPTAKE_RI_EXTMTB);
		cell.setExtMtb(std::max<double>(cell.getExtMtb() - _PARAM(PARAM_MAC_NR_UPTAKE_RI_EXTMTB), 0.0));

		_nextState = MAC_INFECTED;
	}

	// macrophage may become activated
	if (_stat1 && _NFkB)
	{
		_intMtb = 0;
		_activationTime = time;
		_deathTime = time + _PARAM(PARAM_MAC_A_AGE);
		_nextState = MAC_ACTIVE;
	}
}

void Mac::handleInfected(const int time, GrGrid& grid, GrStat&)
{
	GridCell& cell = grid(_row, _col);

	// intracellular bacteria reproduce
	_intMtb *= _PARAM(PARAM_INTMTB_GROWTH_RATE);

	// uptake extracellular bacteria
	if (cell.getExtMtb() > 0 &&
		g_Rand.getReal() < (_PARAM(PARAM_MAC_THRESHOLD_BECOME_CI_INTMTB) - _intMtb) / 100)
	{
		double dExtMtb = std::min(cell.getExtMtb(), _PARAM(PARAM_MAC_NR_UPTAKE_RI_EXTMTB));
		cell.incExtMtb(-dExtMtb);
		_intMtb += dExtMtb;
	}

	if (_intMtb >= _PARAM(PARAM_MAC_THRESHOLD_BECOME_CI_INTMTB))
	{
		_nextState = MAC_CINFECTED;
	}
	else
	{
		_stat1 |= g_Rand.getReal() <
			getCountTgam(TGAM_ACTIVE, grid) * _PARAM(PARAM_MAC_PROB_STAT1_TGAM);

		if (_stat1 && _NFkB)
		{
			_intMtb = 0;
			_activationTime = time;
			_deathTime = time + _PARAM(PARAM_MAC_A_AGE);
			_nextState = MAC_ACTIVE;
		}
		else
		{
			_nextState = MAC_INFECTED;
		}
	}
}

void Mac::handleChronicallyInfected(const int, GrGrid& grid, GrStat&)
{
	GridCell& cell = grid(_row, _col);

	// intracellular bacteria reproduce
	_intMtb *= _PARAM(PARAM_INTMTB_GROWTH_RATE);

	if (_intMtb >= _PARAM(PARAM_MAC_THRESHOLD_BURST_CI_INTMTB))
	{
		// burst, all intracellular bacteria disperse to the Moore neighborhood
		const double dExtMtb = _intMtb / 9.0;

		for (int i = -1; i <= 1; i++)
			for (int j = -1; j <= 1; j++)
				grid(MOD_ROW(_row + i), MOD_COL(_col + j)).incExtMtb(dExtMtb);

		_intMtb = 0;

		// increment number of killings
		cell.incNrKillings();

		_nextState = MAC_DEAD;
	}
	else
	{
		_nextState = MAC_CINFECTED;
	}
}

void Mac::handleActivated(const int, GrGrid& grid, GrStat&)
{
	GridCell& cell = grid(_row, _col);

	// kill extracellular bacteria in the compartment the macrophage resides
	cell.incExtMtb(-1 * std::min(cell.getExtMtb(), _PARAM(PARAM_MAC_NR_UPTAKE_A_EXTMTB)));

	_nextState = MAC_ACTIVE;
}

void Mac::deactivate(const int time)
{
	switch (_state)
	{
	case MAC_RESTING:
	case MAC_INFECTED:
	case MAC_ACTIVE:
		_stat1 = false;
		_deactivationTime = time;
		break;
	default:
		;
	}
}

void Mac::print() const
{
	std::cout << "Mac - ";

	std::cout << _birthTime << " - ";

	switch (_state)
	{
	case MAC_DEAD:
		std::cout << "dead, ";
		break;
	case MAC_RESTING:
		std::cout << "resting, ";
		break;
	case MAC_INFECTED:
		std::cout << "infected, ";
		break;
	case MAC_CINFECTED:
		std::cout << "chronically infected, ";
		break;
	case MAC_ACTIVE:
		std::cout << "active, ";
		break;
	}

	if (_NFkB)
	{
		std::cout << "NFkB, ";
	}

	if (_stat1)
	{
		std::cout << "stat1, ";
	}

	std::cout << _intMtb << std::endl;
}

void Mac::serialize(std::ostream& out) const
{
	assert(out.good());

	Agent::serialize(out);

	int intVal = (int) _state;
	out << intVal << std::endl;

	intVal = (int) _nextState;
	out << intVal << std::endl;

	out << _intMtb << std::endl;
	out << _NFkB << std::endl;
	out << _stat1 << std::endl;
	out << _activationTime << std::endl;
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

void Mac::deserialize(std::istream& in)
{
	assert(in.good());

	int intVal;

	Agent::deserialize(in);

	in >> intVal;
	_state = (MacState) intVal;

	in >> intVal;
	_nextState = (MacState) intVal;

	in >> _intMtb;
	in >> _NFkB;
	in >> _stat1;
	in >> _activationTime;
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

int Mac::getCountTgam(TgamState state, const GrGrid& grid) const
{
	int count = 0;

	// count the number of Tgam cells in the Moore neighborhood
	// that are in the specified state
	for (int i = -1; i <= 1; i++)
		for (int j = -1; j <= 1; j++)
			count += grid(MOD_ROW(_row + i), MOD_COL(_col + j)).hasTgam(state);

	return count;
}

double Mac::getExtMtbInMoore(const GrGrid& grid) const
{
	double mtb = 0;

	// count the number of Tgam cells in the Moore neighborhood
	// that are in the specified state
	for (int i = -1; i <= 1; i++)
		for (int j = -1; j <= 1; j++)
			mtb += grid(MOD_ROW(_row + i), MOD_COL(_col + j)).getExtMtb();

	return mtb;
}
