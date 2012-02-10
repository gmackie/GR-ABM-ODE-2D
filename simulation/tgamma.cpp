/*
 * tgamma.cpp
 *
 *  Created on: 03-nov-2009
 *      Author: M. El-Kebir
 */

#include "tgamma.h"
#include "grgrid.h"
#include "serialization.h"

const std::string Tgam::_ClassName = "Tgam";


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
    , _mTNFRNA(-1.0)
	, _vTNFR1(-1.0)
	, _vTNFR2(-1.0)
	, _kSynth(-1.0)
	, _kTACE(-1.0)
    , _kmRNA(-1.0)

    // IL10 components
    
    , _surfIL10R(-1.0)
    , _vIL10R(-1.0)
    , _surfBoundIL10R(-1.0)
    , _kISynth(-1.0)

    , _nAntigenStim(-1)
    , _nDownRegulated(-1)
    , _nICOS(-1)
{
}

Tgam::Tgam(int birthtime, int row, int col, TgamState state)
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
    , _mTNFRNA(0.0)
	, _vTNFR1(_surfTNFR1 * _PARAM(PARAM_GR_K_T1))
	, _vTNFR2(_surfTNFR2 * _PARAM(PARAM_GR_K_T2))
	, _kSynth(_PARAM(PARAM_GR_K_SYNTH_TCELL))
	, _kTACE(_PARAM(PARAM_GR_K_TACE_TCELL))
    , _kmRNA(0.0)

    // IL10 components
    
    , _surfIL10R(g_Rand.getReal(_PARAM(PARAM_GR_I_IL10R_TCELL) - _PARAM(PARAM_GR_STD_IL10R_TCELL), _PARAM(PARAM_GR_I_IL10R_TCELL) + _PARAM(PARAM_GR_STD_IL10R_TCELL)))
    , _vIL10R(_surfIL10R * _PARAM(PARAM_GR_I_K_T))
    , _surfBoundIL10R(0.0)
    , _kISynth(0.0)

    , _nAntigenStim(0)
    , _nDownRegulated(0)
    , _nICOS(0)
{
}

Tgam::~Tgam()
{
}

void Tgam::move(GrGrid& grid)
{
	Tcell::moveTcell(grid, true, true, true);
}

void Tgam::secrete(GrGrid& grid, bool tnfrDynamics, bool, bool tnfDepletion, bool il10rDynamics, bool il10Depletion)
{
	if (_deactivationTime != -1)
	{
		_kSynth = 0;
        _kmRNA = 0;
        _kISynth = 0;
		return;
	}
	
	GridCell& cell = grid(_row, _col);
	_kSynth = _PARAM(PARAM_GR_K_SYNTH_TCELL);
    _kmRNA = _PARAM(PARAM_GR_K_RNA_TCELL);
    
    
    if (_state == TGAM_ACTIVE_DOUBLE) {
        _kISynth = _PARAM(PARAM_GR_I_K_SYNTH_TCELL);
    }
    else
    {
        _kISynth = 0;
    }
    
    double Nav = 6.02e23; // Avogadro Number
    double vol = 8.0e-12; // volume of a cell in liter
    double MW_IL10 = 18600; // molecular weight of IL10 in g/mol
    
    // Need a statement checking to see if cell is in double positive state
    // Then set _kISynth = T cell value
    
    
	if (!tnfrDynamics && !tnfDepletion)
    {    
        double il10 = log(((cell.getIL10() * MW_IL10 * 1e6)/(Nav * vol))); // converting il10 concentration to log(ng/mL) for use in dose dependence
        double tnfMOD = (1.0/(1.0 + exp((il10 + _PARAM(PARAM_GR_LINK_LOG_ALPHA))/_PARAM(PARAM_GR_LINK_LOG_BETA)))); // calculate the fraction of inhibition
		cell.incTNF(tnfMOD * _PARAM(PARAM_TGAM_SEC_RATE_TNF));
    }
}

void Tgam::computeNextState(const int time, GrGrid& grid, GrStat& stats, bool tnfrDynamics, bool, bool)
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
		stats.incTcellApoptosisTNF();
		_nextState = TGAM_DEAD;
	}
	else if (!tnfrDynamics && tnfBoundFraction > _PARAM(PARAM_GR_THRESHOLD_APOPTOSIS_TNF) &&
			 g_Rand.getReal() < 1 - pow(2.7183, -_PARAM(PARAM_GR_K_APOPTOSIS) * (tnfBoundFraction - _PARAM(PARAM_GR_THRESHOLD_APOPTOSIS_TNF))))
	{
		// TNF induced apoptosis
		stats.incTcellApoptosisTNF();
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
        case TGAM_ACTIVE_DOUBLE:
            handleActiveDouble(time, grid, stats);
            break;
                
    default: throw std::runtime_error("Unknown Tgam state");
		}
	}
}

void Tgam::handleActive(const int, GrGrid& grid, GrStat& stats)
{
	GridCell& cell = grid(_row, _col);
	
    double ProbSum = 0.0; // Initialize the probability sum for transition to TGAM_ACTIVE_DOUBLE

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
        
        
        // If mac does not die then check for Tgam10 conditions - Macs
        if (pMac && (pMac->getState() == MAC_INFECTED || pMac->getState() == MAC_CINFECTED))
        {
            
            if (pMac && pMac->getStat1())
                {
                    _nICOS ++;
                    cout << "ICOS" << std::endl;
                }
            
            if (g_Rand.getReal() < _PARAM(PARAM_TGAM_PROB_ANTIGEN_PRESENTATION))
                {
                    _nAntigenStim ++;
                    cout << "Antigen Stim" << std::endl;
                }
            
        }
        if (pMac && ((pMac->getState() == MAC_RESTING && pMac->getNFkB()) || (pMac->getState() == MAC_RESTING && pMac->getStat1()) || (pMac->getState() == MAC_ACTIVE)) && cell.getExtMtb() > _PARAM(PARAM_TGAM_THRESHOLD_EXT_MTB))
        {
            _nAntigenStim ++;
            cout << "Antigen Stim" << std::endl;
        }
        if (pMac && ((pMac->getState() == MAC_RESTING && pMac->getStat1()) || (pMac->getState() == MAC_ACTIVE && pMac->getStat1())))
        {
            _nICOS ++;
            cout << "ICOS" << std::endl;
        }
        
	}
    
    if (_nAntigenStim > 1) 
    {
        
        double AgProb = 0.5 * (1 - pow(2.7183, (-_PARAM(PARAM_TGAM_PROB_AGSTIM)*(_nAntigenStim - 1.0))));
        ProbSum += AgProb;
        
    }
    
    
    if (_nDownRegulated > 0) 
    {
        
        double TGFProb = 0.25 * (1 - pow(2.7183, (-_PARAM(PARAM_TGAM_PROB_AGSTIM)*(_nDownRegulated))));
        ProbSum += TGFProb;
        
    }
    
    if (_nICOS > 0) 
    {
        
        double ICOSProb = 0.25 * (1 - pow(2.7183, (-_PARAM(PARAM_TGAM_PROB_AGSTIM)*(_nICOS))));
        ProbSum += ICOSProb;
        
    }
    
    cout << "ProbabilitySum:   " << ProbSum << std::endl;
    
    
    if (g_Rand.getReal() < ProbSum) {
        
        _nextState = TGAM_ACTIVE_DOUBLE;
    }
    else
    {
        _nextState = TGAM_ACTIVE;
    }
    
    if (ProbSum < 0 || ProbSum > 1) {
        std::cout << "Error:Probability of Transition to IL10 Producing Tgamma is  " << ProbSum << std::endl;
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

void Tgam::handleActiveDouble(const int time, GrGrid& grid, GrStat& stats)
{
    // Carries out same action as TGAM_ACTIVE but secretes IL10
    // Distinct state so it is easier to track
    
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
        
		_nextState = TGAM_ACTIVE_DOUBLE;
        
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

void Tgam::deactivate(const int time)
{
    if (_state == TGAM_ACTIVE || _state == TGAM_ACTIVE_DOUBLE) {
        _nDownRegulated ++;
        cout << _nDownRegulated << std::endl;
        cout << "TGF-B" << std::endl;
    }
	_nextState = _state = TGAM_DOWN_REGULATED;
	_deactivationTime = time;

}

void Tgam::updateState()
{
	_state = _nextState;
}

void Tgam::solveTNF(GrGrid& grid, double dt)
{
    GridCell& cell = grid(_row, _col);
    
	double koff1 = _PARAM(PARAM_GR_K_ON1) * _PARAM(PARAM_GR_KD1);
	double koff2 = _PARAM(PARAM_GR_K_ON2) * _PARAM(PARAM_GR_KD2);
	double density = 1.25e11; // used for conversion of conc. unit (M -> #/cell) based on cell and microcompartment volumes 
	double Nav = 6.02e23; // Avogadro Number
	double vol = 8.0e-12; // volume of a cell in liter
	
	double tnf = cell.getTNF() / (Nav * vol);
	double shedtnfr2 = cell.getShedTNFR2() / (Nav * vol);
    double il10 = cell.getIL10() /(Nav * vol);
	
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
	
    double eqsurfBoundIL10R;
    double IkmRNA;
    
    // modulate TNF parameters based on equilibrium IL10 receptor equations since the il10 molecular equations are not active
    eqsurfBoundIL10R = (il10 * _surfIL10R) / (_PARAM(PARAM_GR_I_KD) + il10);
    if (_kmRNA > 0) {
        IkmRNA = _kmRNA * ((_kSynth/_kmRNA)+ ((1.0 - (_kSynth/_kmRNA))/(1.0 + pow(2.7183, ((eqsurfBoundIL10R - _PARAM(PARAM_GR_LINK_RNA_GAMMA))/_PARAM(PARAM_GR_LINK_RNA_DELTA))))));
    }
    else
    {
        IkmRNA = 0.0;
    }
    // end of equilibrium calculations
    
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
	
	cell.setTNF(Nav * vol * tnf);
	cell.setShedTNFR2(Nav * vol * shedtnfr2);
	if (_mTNF < 0 || _surfTNFR1 < 0 || _surfBoundTNFR1 < 0 || _surfTNFR2 < 0 || _surfBoundTNFR2 < 0 || _mTNFRNA < 0)
		std::cout << "Error: Negative Value of Species in TNF/TNFR dynamics" << std::endl;
    
    //cout << "Debug: Running TNF dynamics" << std::endlcies in TNF/TNFR dynamics" << std::endl;

}

void Tgam::solveTNFandIL10(GrGrid& grid, double dt)
{
    GridCell& cell = grid(_row, _col);
	
	double koff1 = _PARAM(PARAM_GR_K_ON1) * _PARAM(PARAM_GR_KD1);
	double koff2 = _PARAM(PARAM_GR_K_ON2) * _PARAM(PARAM_GR_KD2);
	double density = 1.25e11; // used for conversion of conc. unit (M -> #/cell) based on cell and microcompartment volumes 
	double Nav = 6.02e23; // Avogadro Number
	double vol = 8.0e-12; // volume of a cell in liter
	
	double tnf = cell.getTNF() / (Nav * vol);
	double shedtnfr2 = cell.getShedTNFR2() / (Nav * vol);
    double il10 = cell.getIL10() / (Nav * vol);
	
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
	
	cell.setTNF(Nav * vol * tnf);
	cell.setShedTNFR2(Nav * vol * shedtnfr2);
    cell.setIL10(Nav * vol * il10);
	
    
	if (_mTNF < 0 || _surfTNFR1 < 0 || _surfBoundTNFR1 < 0 || _surfTNFR2 < 0 || _surfBoundTNFR2 < 0 || _mTNFRNA < 0)
		std::cout << "Error: Negative Value of Species in TNF/TNFR dynamics" << std::endl;
    
    if (_surfIL10R < 0 || _surfBoundIL10R < 0)
        std::cout << "Error: Negative value of species in IL10/IL10R dynamics" << std::endl;
    
    //cout << "Debug: Running TNF and IL10 dynamics" << std::endl;    
}

void Tgam::solveIL10(GrGrid& grid, double dt)
{
    GridCell& cell = grid(_row, _col);
    
    double density = 1.25e11; // used for conversion of conc. unit (M -> #/cell) based on cell and microcompartment volumes 
	double Nav = 6.02e23; // Avogadro Number
	double vol = 8.0e-12; // volume of a cell in liter
    
    double il10 = cell.getIL10() / (Nav * vol);
    
    double dsIL10;
	double dsurfIL10R;
	double dsurfBoundIL10R;
    
    // IL10 differential equations
    dsIL10 = (density/Nav) * _kISynth + ((density/Nav) * (_PARAM(PARAM_GR_I_K_OFF) * _surfBoundIL10R - _PARAM(PARAM_GR_I_K_ON) * _surfIL10R * il10)) * dt;
	dsurfIL10R = (_vIL10R - _PARAM(PARAM_GR_I_K_ON) * _surfIL10R * il10 + _PARAM(PARAM_GR_I_K_OFF) * _surfBoundIL10R - _PARAM(PARAM_GR_I_K_T) * _surfIL10R) * dt;
	dsurfBoundIL10R = (_PARAM(PARAM_GR_I_K_ON) * _surfIL10R * il10 - _PARAM(PARAM_GR_I_K_OFF) * _surfBoundIL10R - _PARAM(PARAM_GR_I_K_INT) * _surfBoundIL10R) * dt;
    // end of IL10 differential equations
    
    // update il10 variables
    _surfIL10R += dsurfIL10R;
    _surfBoundIL10R += dsurfBoundIL10R;
    il10 += dsIL10;
    
    cell.setIL10(Nav * vol * il10);
    
    if (_surfIL10R < 0 || _surfBoundIL10R < 0)
        std::cout << "Error: Negative value of species in IL10/IL10R dynamics" << std::endl;
    
}

void Tgam::solveDegradation(GrGrid& grid, double dt, bool tnfrDynamics, bool il10rDynamics)
{
    GridCell& cell = grid(_row, _col);
    
    double Nav = 6.02e23; // Avagadros #
    double vol = 8.0e-12; // volume of a cell in L
    
    if (!tnfrDynamics) {
        
        // simulate the effect of TNF internalization by cells in the form of degradation. Only for TNF
        double dtnf;
        double tnf = cell.getTNF();
        dtnf = -_PARAM(PARAM_GR_K_INT1) * (tnf / (tnf + _PARAM(PARAM_GR_KD1) * Nav * vol)) * _PARAM(PARAM_GR_MEAN_TNFR1_TCELL) * dt * 0.4;
        tnf += dtnf;  
        
        cell.setTNF(tnf);
    }
    
    if (!il10rDynamics) {
        
        double dil10;
        double il10 = cell.getIL10();
        
        // simulate the effect of IL10 internalization in the form of degradation. Only for IL10
        dil10 = -_PARAM(PARAM_GR_I_K_INT) * (il10 / (il10 + _PARAM(PARAM_GR_I_KD) * Nav * vol)) * _PARAM(PARAM_GR_I_IL10R_TCELL) * dt * _PARAM(PARAM_GR_I_MOD);
        il10 += dil10;  
        
        cell.setIL10(il10);
        
    }
    
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
  default: throw std::runtime_error("Unknown Tcyt state");
	}
}

void Tgam::serialize(std::ostream& out) const
{
	assert(out.good());

	Serialization::writeHeader(out, Tgam::_ClassName);

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
    out << _mTNFRNA << std::endl;
	out << _vTNFR1 << std::endl;
	out << _vTNFR2 << std::endl;
	out << _kSynth << std::endl;
	out << _kTACE << std::endl;
    out << _kmRNA << std::endl;
    
    out << _surfIL10R << std::endl;
    out << _vIL10R << std::endl;
    out << _surfBoundIL10R << std::endl;
    out << _kISynth << std::endl;
    
    out << _nAntigenStim << std::endl;
    out << _nDownRegulated << std::endl;
    out << _nICOS << std::endl;
	
	Serialization::writeFooter(out, Tgam::_ClassName);
}

void Tgam::deserialize(std::istream& in)
{
	assert(in.good());

	if (!Serialization::readHeader(in, Tgam::_ClassName))
	{
		exit(1);
	}

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
    in >> _mTNFRNA;
	in >> _vTNFR1;
	in >> _vTNFR2;
	in >> _kSynth;
	in >> _kTACE;
    in >> _kmRNA;
    
    in >> _surfIL10R;
    in >> _vIL10R;
    in >> _surfBoundIL10R;
    in >> _kISynth;
    
    in >> _nAntigenStim;
    in >> _nDownRegulated;
    in >> _nICOS;

	if (!Serialization::readFooter(in, Tgam::_ClassName))
	{
		exit(1);
	}
}
