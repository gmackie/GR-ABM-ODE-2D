/*
 * tgamma.cpp
 *
 *  Created on: 03-nov-2009
 *      Author: M. El-Kebir
 */

#include "tgamma.h"
#include "macrophage.h"
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
    , _transitionTime(-1)
    , _nAntigenStim(-1)
    , _nDownRegulated(-1)
    , _nICOS(-1)
{
}

Tgam::Tgam(int birthtime, int row, int col, TgamState state)

	: Tcell(birthtime, row, col, _PARAM(PARAM_GR_K_SYNTH_TCELL))
	, _state(state)
	, _nextState(state)
	, _deactivationTime(-1)
    , _transitionTime(-1)

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
	
	_kSynth = _PARAM(PARAM_GR_K_SYNTH_TCELL);
    _kmRNA = _PARAM(PARAM_GR_K_RNA_TCELL);
    
    
    if (_state == TGAM_ACTIVE_DOUBLE) {
        _kISynth = 0.5 * _PARAM(PARAM_GR_I_K_SYNTH_TCELL);
    }
    
    else if (_state == TGAM_INDUCED_REG)
    {
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
        double il10 = log(((grid.il10(_pos) * MW_IL10 * 1e6)/(Nav * vol))); // converting il10 concentration to log(ng/mL) for use in dose dependence
        double tnfMOD = (1.0/(1.0 + exp((il10 + _PARAM(PARAM_GR_LINK_LOG_ALPHA))/_PARAM(PARAM_GR_LINK_LOG_BETA)))); // calculate the fraction of inhibition
		grid.incTNF(_pos, (tnfMOD * _PARAM(PARAM_TGAM_SEC_RATE_TNF)));
    }
}

void Tgam::computeNextState(const int time, GrGrid& grid, GrStat& stats, bool tnfrDynamics, bool, bool, bool tgammatransition)
{
    
	double tnfBoundFraction = grid.TNF(_pos) / (grid.TNF(_pos) + _PARAM(PARAM_GR_KD1) * 48.16e11);

    
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
			handleActive(time, grid, stats, tgammatransition);
			break;
		case TGAM_DOWN_REGULATED:
			handleDownRegulated(time, grid, stats);
			break;
        case TGAM_ACTIVE_DOUBLE:
            handleActiveDouble(time, grid, stats);
            break;
        case TGAM_INDUCED_REG:
            handleInducedReg(time, grid, stats);
            break;
                
    default: throw std::runtime_error("Unknown Tgam state");
		}
	}
}

void Tgam::handleActive(const int time, GrGrid& grid, GrStat& stats, bool tgammatransition)
{
    double ProbSum = 0.0; // Initialize the probability sum for transition to TGAM_ACTIVE_DOUBLE

    if(tgammatransition)
	{
		if (grid.hasAgentType(MAC, _pos))
		{
			// get the macrophage
			Mac* pMac = dynamic_cast<Mac*>(grid.agent(_pos, 0));
			if (!pMac) pMac = dynamic_cast<Mac*>(grid.agent(_pos, 1));

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

				grid.incKillings(_pos);
			}


            
            // If mac does not die then check for Tgam10 conditions - Macs
            if (pMac && (pMac->getState() == MAC_INFECTED || pMac->getState() == MAC_CINFECTED))
            {
                
                if (pMac && pMac->getICOS())
                {
                    _nICOS ++;
//                    cout << "ICOS" << std::endl;
                }
                
                if (g_Rand.getReal() < _PARAM(PARAM_TGAM_PROB_ANTIGEN_PRESENTATION))
                {
                    _nAntigenStim ++;
//                    cout << "Antigen Stim" << std::endl;
                }
                
            }
        	if (pMac && ((pMac->getState() == MAC_RESTING && pMac->getNFkB()) || (pMac->getState() == MAC_RESTING && pMac->getICOS()) || (pMac->getState() == MAC_ACTIVE)) && grid.extMTB(_pos) > _PARAM(PARAM_TGAM_THRESHOLD_EXT_MTB))
        	{
         	   _nAntigenStim ++;
//         	   cout << "Antigen Stim" << std::endl;
       		 }
        	if (pMac && ((pMac->getState() == MAC_RESTING && pMac->getICOS()) || (pMac->getState() == MAC_ACTIVE && pMac->getICOS())))
       	 {
	            _nICOS ++;
// 	           cout << "ICOS" << std::endl;
	        }
            
        }
        
        if (_nAntigenStim > 1) 
        {
            
            double AgProb = (1.0/3.0) * (1 - pow(2.7183, (-_PARAM(PARAM_TGAM_PROB_AGSTIM)*(_nAntigenStim - 1.0))));
            ProbSum += AgProb;
            
        }
        
        
        if (_nDownRegulated > 0) 
        {
            
            double TGFProb = (1.0/3.0) * (1 - pow(2.7183, (-_PARAM(PARAM_TGAM_PROB_TGFB)*(_nDownRegulated))));
            ProbSum += TGFProb;
            
        }
        
        if (_nICOS > 0) 
        {
            
            double ICOSProb = (1.0/3.0) * (1 - pow(2.7183, (-_PARAM(PARAM_TGAM_PROB_ICOS)*(_nICOS))));
            
            ProbSum += ICOSProb;
            
        }
        
//        cout << "ProbabilitySum:   " << ProbSum << std::endl;
        
        
        if (g_Rand.getReal() < ProbSum) {
            
            _nextState = TGAM_ACTIVE_DOUBLE;
            _transitionTime = time;
        }
        else
        {
            _nextState = TGAM_ACTIVE;
        }
        
        if (ProbSum < 0 || ProbSum > 1) {
            std::cout << "Error: Probability of Transition to IL10 Producing Tgamma is  " << ProbSum << std::endl;
        }

    }
    
    else
    {
        if (grid.hasAgentType(MAC, _pos))
        {
            // get the macrophage
            Mac* pMac = dynamic_cast<Mac*>(grid.agent(_pos, 0));
            if (!pMac) pMac = dynamic_cast<Mac*>(grid.agent(_pos, 1));
            
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
                
                grid.incKillings(_pos);
            }
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

void Tgam::handleActiveDouble(const int time, GrGrid& grid, GrStat& stats)
{
    // Carries out same action as TGAM_ACTIVE but secretes half rate of IL10
    // Distinct state so it is easier to track
    
    if (time - _transitionTime >= _PARAM(PARAM_TGAM_TIMESPAN_DOUBLE)) {
        _nextState = TGAM_INDUCED_REG;
    }
    else
    {
        _nextState = TGAM_ACTIVE_DOUBLE;
    }

	if (grid.hasAgentType(MAC, _pos))
	{
		// get the macrophage
		Mac* pMac = dynamic_cast<Mac*>(grid.agent(_pos, 0));
		if (!pMac) pMac = dynamic_cast<Mac*>(grid.agent(_pos, 1));
        
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
            
			grid.incKillings(_pos);
		}
	}
}

void Tgam::handleInducedReg(const int time, GrGrid& grid, GrStat& stats)
{
    _nextState = TGAM_INDUCED_REG;
}

void Tgam::deactivate(const int time)
{
    if (_state == TGAM_ACTIVE || _state == TGAM_ACTIVE_DOUBLE) {
        _nDownRegulated ++;
//        cout << _nDownRegulated << std::endl;
//        cout << "TGF-B" << std::endl;
    }
    
    if (_state == TGAM_INDUCED_REG) {
        _nextState = TGAM_INDUCED_REG;
    }
    else
    {
	_nextState = _state = TGAM_DOWN_REGULATED;
	_deactivationTime = time;
    }
}

void Tgam::updateState()
{
	_state = _nextState;
}

void Tgam::solveDegradation(GrGrid& grid, double dt, bool tnfrDynamics, bool il10rDynamics)
{
    double Nav = 6.02e23; // Avagadros #
    double vol = 8.0e-12; // volume of a cell in L
    
    if (!tnfrDynamics) {
        
        // simulate the effect of TNF internalization by cells in the form of degradation. Only for TNF
        double dtnf;
        double tnf = grid.TNF(_pos);
        dtnf = -_PARAM(PARAM_GR_K_INT1) * (tnf / (tnf + _PARAM(PARAM_GR_KD1) * Nav * vol)) * _PARAM(PARAM_GR_MEAN_TNFR1_TCELL) * dt * 0.4;
        grid.incTNF(_pos, dtnf);  
    }
    
    if (!il10rDynamics) {
        
        double dil10;
        double il10 = grid.il10(_pos);
        
        // simulate the effect of IL10 internalization in the form of degradation. Only for IL10
        dil10 = -_PARAM(PARAM_GR_I_K_INT) * (il10 / (il10 + _PARAM(PARAM_GR_I_KD) * Nav * vol)) * _PARAM(PARAM_GR_I_IL10R_TCELL) * dt * _PARAM(PARAM_GR_I_MOD);
        grid.incil10(_pos, dil10);  
        
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
    out << _transitionTime << std::endl;

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
    in >> _transitionTime;
    
    in >> _nAntigenStim;
    in >> _nDownRegulated;
    in >> _nICOS;

	if (!Serialization::readFooter(in, Tgam::_ClassName))
	{
		exit(1);
	}
}
