/*
 * tcell.cpp
 *
 *  Created on: 13-nov-2009
 *      Author: M. El-Kebir
 */


#include "tcell.h"
#include "grgrid.h"

using namespace std;

// Needed for deserializing the model state.
// Avoids the calls to the random number generator in the normal constructor, allowing the random number generator
// to remain in synch after deserialization.
Tcell::Tcell()
	: Agent()
{
}

Tcell::Tcell(int birthtime, int row, int col, Scalar kSynth)
	: Agent(birthtime, birthtime + _PARAM(PARAM_TCELL_AGE), row, col
			//TNFR Components
			, (Scalar) _PARAM(PARAM_GR_MEAN_TNFR1_TCELL)
			, (Scalar) _PARAM(PARAM_GR_STD_TNFR1_TCELL)
			, (Scalar) _PARAM(PARAM_GR_MEAN_TNFR2_TCELL)
			, (Scalar) _PARAM(PARAM_GR_STD_TNFR2_TCELL)
			, kSynth
			, (Scalar) _PARAM(PARAM_GR_K_TACE_TCELL)

			// IL10 components
			, (Scalar) _PARAM(PARAM_GR_I_IL10R_TCELL)
			, (Scalar) _PARAM(PARAM_GR_STD_IL10R_TCELL)
			)
{
}

Tcell::~Tcell()
{
}

void Tcell::solveTNFandIL10(GrGrid& grid, GrStat&, double dt, double)
{
	double koff1 = _PARAM(PARAM_GR_K_ON1) * _PARAM(PARAM_GR_KD1);
	double koff2 = _PARAM(PARAM_GR_K_ON2) * _PARAM(PARAM_GR_KD2);
	double density = 1.25e11; // used for conversion of conc. unit (M -> #/cell) based on cell and microcompartment volumes
	double Nav = 6.02e23; // Avogadro Number
	double vol = 8.0e-12; // volume of a cell in liter

	double tnf = grid.TNF(_pos) / (Nav * vol);
	double shedtnfr2 = grid.shedTNFR2(_pos) / (Nav * vol);
    double il10 = grid.il10(_pos) / (Nav * vol);

    double dmTNFRNA;
	double dmTNF;
	double dsurfTNFR1;
	double dsurfTNFR2;
	double dsurfBoundTNFR1;
	double dsurfBoundTNFR2;
	double dintBoundTNFR1;
	double dintBoundTNFR2;
	double dsTNF;
	double dshedTNFR2;

    double dsIL10;
	double dsurfIL10R;
	double dsurfBoundIL10R;

    double IkmRNA;

    // solving for TNF parameters that depend on IL10

    if (_kmRNA > 0) {
        IkmRNA = _kmRNA * ((_kSynth/_kmRNA)+ ((1.0 - (_kSynth/_kmRNA))/(1.0 + pow(2.7183, ((_surfBoundIL10R - _PARAM(PARAM_GR_LINK_RNA_GAMMA))/_PARAM(PARAM_GR_LINK_RNA_DELTA))))));
    }
    else
    {
        IkmRNA = 0.0;
    }

    // end of TNF and IL10 linking

	// TNF differential equations
	dmTNFRNA = (IkmRNA - _PARAM(PARAM_GR_K_TRANS) * _mTNFRNA) * dt;
	dmTNF = (_PARAM(PARAM_GR_K_TRANS) * _mTNFRNA - _kTACE * _mTNF) * dt;
	dsurfTNFR1 = (_vTNFR1 - _PARAM(PARAM_GR_K_ON1) * tnf * _surfTNFR1 + koff1 * _surfBoundTNFR1 - _PARAM(PARAM_GR_K_T1) * _surfTNFR1 + _PARAM(PARAM_GR_K_REC1) * _intBoundTNFR1) * dt;
	dsurfTNFR2 = (_vTNFR2 - _PARAM(PARAM_GR_K_ON2) * tnf * _surfTNFR2 + koff2 * _surfBoundTNFR2 - _PARAM(PARAM_GR_K_T2) * _surfTNFR2 + _PARAM(PARAM_GR_K_REC2) * _intBoundTNFR2) * dt;
	dsurfBoundTNFR1 = (_PARAM(PARAM_GR_K_ON1) * tnf * _surfTNFR1 - koff1 * _surfBoundTNFR1 - _PARAM(PARAM_GR_K_INT1) * _surfBoundTNFR1) * dt;
	dsurfBoundTNFR2 = (_PARAM(PARAM_GR_K_ON2) * tnf * _surfTNFR2 - koff2 * _surfBoundTNFR2 - _PARAM(PARAM_GR_K_INT2) * _surfBoundTNFR2 - _PARAM(PARAM_GR_K_SHED) * _surfBoundTNFR2) * dt;
	dintBoundTNFR1 = (_PARAM(PARAM_GR_K_INT1) * _surfBoundTNFR1 - _PARAM(PARAM_GR_K_DEG1) * _intBoundTNFR1 - _PARAM(PARAM_GR_K_REC1) * _intBoundTNFR1) * dt;
	dintBoundTNFR2 = (_PARAM(PARAM_GR_K_INT2) * _surfBoundTNFR2 - _PARAM(PARAM_GR_K_DEG2) * _intBoundTNFR2 - _PARAM(PARAM_GR_K_REC2) * _intBoundTNFR2) * dt;
	dsTNF = ((density/Nav) * (_kTACE * _mTNF - _PARAM(PARAM_GR_K_ON1) * tnf * _surfTNFR1 + koff1 * _surfBoundTNFR1 - _PARAM(PARAM_GR_K_ON2) * tnf * _surfTNFR2 + koff2 * _surfBoundTNFR2)) * dt;
	dshedTNFR2 = ((density/Nav) * _PARAM(PARAM_GR_K_SHED) * _surfBoundTNFR2) * dt;
	// end of TNF differential equations

    // IL10 differential equations
    dsIL10 = (density/Nav) * _kISynth + ((density/Nav) * (_PARAM(PARAM_GR_I_K_OFF) * _surfBoundIL10R - _PARAM(PARAM_GR_I_K_ON) * _surfIL10R * il10)) * dt;
	dsurfIL10R = (_vIL10R - _PARAM(PARAM_GR_I_K_ON) * _surfIL10R * il10 + _PARAM(PARAM_GR_I_K_OFF) * _surfBoundIL10R - _PARAM(PARAM_GR_I_K_T) * _surfIL10R) * dt;
	dsurfBoundIL10R = (_PARAM(PARAM_GR_I_K_ON) * _surfIL10R * il10 - _PARAM(PARAM_GR_I_K_OFF) * _surfBoundIL10R - _PARAM(PARAM_GR_I_K_INT) * _surfBoundIL10R) * dt;
    // end of IL10 differential equations

    // update tnf variables
	_mTNFRNA += dmTNFRNA;
	_mTNF += dmTNF;
	_surfTNFR1 += dsurfTNFR1;
	_surfTNFR2 += dsurfTNFR2;
	_surfBoundTNFR1 += dsurfBoundTNFR1;
	_surfBoundTNFR2 += dsurfBoundTNFR2;
	_intBoundTNFR1 += dintBoundTNFR1;
	_intBoundTNFR2 += dintBoundTNFR2;
	tnf += dsTNF;
	shedtnfr2 += dshedTNFR2;

    // update il10 variables
    _surfIL10R += dsurfIL10R;
    _surfBoundIL10R += dsurfBoundIL10R;
    il10 += dsIL10;

    grid.setTNF(_pos, (Nav * vol * tnf));
    grid.setshedTNFR2(_pos, (Nav * vol * shedtnfr2));
    grid.setil10(_pos, (Nav * vol * il10));


	if (_mTNF < 0 || _surfTNFR1 < 0 || _surfBoundTNFR1 < 0 || _surfTNFR2 < 0 || _surfBoundTNFR2 < 0 || _mTNFRNA < 0)
		std::cout << "Error: Negative Value of Species in TNF/TNFR dynamics" << std::endl;

    if (_surfIL10R < 0 || _surfBoundIL10R < 0)
        std::cout << "Error: Negative value of species in IL10/IL10R dynamics" << std::endl;

//    cout << "Debug: Running T cell TNF and IL10 dynamics" << std::endl;
}

void Tcell::moveTcell(GrGrid& grid, bool ccl2, bool ccl5, bool cxcl9)
{
	Pos pos  = Agent::moveAgent(grid, ccl2, ccl5, cxcl9, false, _PARAM(PARAM_TCELL_MOVEMENT_BONUSFACTOR));
	
	// Check whether newCell is not caseated and contains empty slots
  size_t n = grid.getNumberOfAgents(pos);
	if (n != 2 && !grid.isCaseated(pos))
	{
		// Move with p = 1, if newCell is empty
		// Move with p = _PROB_MOVE_TCELL_TO_MAC, if newCell contains a macrophage
		// Move with p = _PROB_MOVE_TCELL_TO_TCELL, if newCell contains a T cell
		if ((n == 0) ||
			(grid.hasAgentType(MAC, pos) && g_Rand.getReal() < _PARAM(PARAM_TCELL_PROB_MOVE_TO_MAC)) ||
			(grid.hasTcell(pos) && g_Rand.getReal() < _PARAM(PARAM_TCELL_PROB_MOVE_TO_TCELL)))
		{
			assert_res(grid.removeAgent(this));
			grid.addAgent(this, pos);
      _pos = pos;
		}
	}
}

void Tcell::serialize(std::ostream& out) const
{
	Agent::serialize(out);
}

void Tcell::deserialize(std::istream& in)
{
	Agent::deserialize(in);
}
