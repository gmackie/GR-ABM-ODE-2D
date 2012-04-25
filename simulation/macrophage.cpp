/*
 * macrophage.cpp
 *
 *  Created on: 05-nov-2009
 *      Author: M. El-Kebir
 */

#include "macrophage.h"
#include "params.h"
#include "grstat.h"
#include "grgrid.h"
#include "serialization.h"

using namespace std;

const std::string Mac::_ClassName = "Mac";

int Mac::_macodeSize = 0;

// Needed for deserializing the model state.
// Avoids the calls to the random number generator in the normal constructor, allowing the random number generator
// to remain in synch after deserialization.
Mac::Mac()
	:Agent()
	, _state(Mac::MAC_DEAD)
	, _nextState(Mac::MAC_DEAD)
	, _intMtb(-1.0)
	, _NFkB(0)
	, _stat1(0)
	, _ICOS(0)
	, _activationTime(-1)
	, _deactivationTime(-1)

{
}

Mac::Mac(int birthtime, int row, int col, Mac::State state, double intMtb, bool NFkB, bool stat1)
	: Agent(birthtime, birthtime + _PARAM(PARAM_MAC_AGE), row, col

			//TNFR Components
			, (Scalar) _PARAM(PARAM_GR_MEAN_TNFR1_MAC)
			, (Scalar) _PARAM(PARAM_GR_STD_TNFR1_MAC)
			, (Scalar) _PARAM(PARAM_GR_MEAN_TNFR2_MAC)
			, (Scalar) _PARAM(PARAM_GR_STD_TNFR2_MAC)
			, 0.0 // kSynth
			, (Scalar) _PARAM(PARAM_GR_K_TACE_MAC)

			// IL10 components
			, (Scalar) _PARAM(PARAM_GR_I_IL10R_MAC)
			, (Scalar) _PARAM(PARAM_GR_STD_IL10R_MAC)
            , _macodeSize
			)
	, _state(state)
	, _nextState(state)
	, _intMtb(intMtb)
	, _NFkB(NFkB)
	, _stat1(stat1)
	, _ICOS(0)
	, _activationTime(-1)
	, _deactivationTime(-1)

{
}

Mac::~Mac()
{
}

void Mac::move(GrGrid& grid)
{
	Pos pos = Agent::moveAgent(grid, true, true, false, true, _PARAM(PARAM_MAC_MOVEMENT_BONUSFACTOR));

	if (grid.addAgent(this, pos))
	{
		assert_res(grid.removeAgent(this));
    _pos = pos;
	}
}

void Mac::secrete(GrGrid& grid, bool tnfrDynamics, bool nfkbDynamics, bool tnfDepletion, bool il10rDynamics, bool il10Depletion, int mdt)
{
	if (_deactivationTime != -1)
	{
		_kSynth = 0;
        _kmRNA = 0;
		_c1rChem = 0;
		_c1rTNF = 0;
		_c1rrChemTNF = 0;
        _kISynth = 0;
		return;
	}
	
	if (nfkbDynamics) // TNF and chemokines are secreted as a function of NFkB dynamics
	{
		if (_state == Mac::MAC_RESTING)
		{
			_c1rChem = _PARAM(PARAM_GR_epsilon2) * _PARAM(PARAM_GR_c1r);
			_c1rTNF = _PARAM(PARAM_GR_epsilon2) * _PARAM(PARAM_GR_c1r);
			_c1rrChemTNF = 0;
            _kISynth = 0.0;
		}
		else if (_state == Mac::MAC_INFECTED)
		{
			_c1rChem = _PARAM(PARAM_GR_c1r);
			_c1rTNF = _PARAM(PARAM_GR_c1r);
			_c1rrChemTNF = _PARAM(PARAM_GR_epsilon1) * _PARAM(PARAM_GR_c1r);
			++grid.nSecretions(_pos);
            _kISynth = _PARAM(PARAM_GR_I_K_SYNTH_MAC_INF);
            
            if (!il10rDynamics && !il10Depletion) {
                grid.incil10(_pos, _PARAM(PARAM_MAC_SEC_RATE_IL10) * mdt);
            }
            
		}
		else if (_state == Mac::MAC_ACTIVE)
		{
			_c1rChem = _PARAM(PARAM_GR_c1r);
			_c1rTNF = _PARAM(PARAM_GR_c1r);
			_c1rrChemTNF = _PARAM(PARAM_GR_c1r);
      ++grid.nSecretions(_pos);
            _kISynth = _PARAM(PARAM_GR_I_K_SYNTH_MAC_ACT);
            
            if (!il10rDynamics && !il10Depletion) {
                grid.incil10(_pos,  0.5 * _PARAM(PARAM_MAC_SEC_RATE_IL10) * mdt);
            }
            
		}
        else if (_state == Mac::MAC_CINFECTED)
        {
            _c1rChem = _PARAM(PARAM_GR_c1r);
			_c1rTNF = _PARAM(PARAM_GR_c1r);
			_c1rrChemTNF = _PARAM(PARAM_GR_c1r);
      ++grid.nSecretions(_pos);
            _kISynth = 2.0 * _PARAM(PARAM_GR_I_K_SYNTH_MAC_INF);
            
            if (!il10rDynamics && !il10Depletion) {
                grid.incil10(_pos, 2.0 * _PARAM(PARAM_MAC_SEC_RATE_IL10) * mdt);
            }
        }
	}
    
	else
	// TNF, IL10, and chemokines are secreted independent of NFkB dynamics
    // IL10 dynamics are independent of the state of _NFkB
	{
		if (_state == Mac::MAC_RESTING) {
            
            if (_NFkB) {
                grid.incCCL2(_pos, (0.5 * _PARAM(PARAM_MAC_SEC_RATE_CCL2) * mdt));
				grid.incCCL5(_pos, (0.5 * _PARAM(PARAM_MAC_SEC_RATE_CCL5) * mdt));
				grid.incCXCL9(_pos,  (0.5 * _PARAM(PARAM_MAC_SEC_RATE_CXCL9) * mdt));
				_kSynth = 0.5 * _PARAM(PARAM_GR_K_SYNTH_MAC);
                _kmRNA = 0.5 * _PARAM(PARAM_GR_K_RNA_MAC);
                _kISynth = 0.0;
                
                if (!tnfrDynamics && !tnfDepletion)
                {
                    
                    double il10 = log(((grid.il10(_pos) * MW_IL10 * 1e6)/(NAV * VOL))); // converting il10 concentration to log(ng/mL) for use in dose dependence
                    //std::cout << "LOG OF IL10: " << il10 << std::endl;
                    double tnfMOD = (1.0/(1.0 + exp((il10 + _PARAM(PARAM_GR_LINK_LOG_ALPHA))/_PARAM(PARAM_GR_LINK_LOG_BETA)))); // calculate the fraction of inhibition
                    //std::cout << "tnfMOD: " << tnfMOD << std::endl;
		    grid.incTNF(_pos, (tnfMOD * 0.5 * _PARAM(PARAM_MAC_SEC_RATE_TNF) * mdt));
                }
                
                ++grid.nSecretions(_pos);
            }
            else
            {
                _kSynth = 0.0;
                _kmRNA = 0.0;
                _kISynth = 0.0;
            }
        }
        else if (_state == Mac::MAC_INFECTED) {
            
            if (_NFkB) {
                grid.incCCL2(_pos, (_PARAM(PARAM_MAC_SEC_RATE_CCL2) * mdt));
                grid.incCCL5(_pos, (_PARAM(PARAM_MAC_SEC_RATE_CCL5) * mdt));
                grid.incCXCL9(_pos,  (_PARAM(PARAM_MAC_SEC_RATE_CXCL9) * mdt));
                _kSynth = _PARAM(PARAM_GR_K_SYNTH_MAC);
                _kmRNA = _PARAM(PARAM_GR_K_RNA_MAC);
                _kISynth = _PARAM(PARAM_GR_I_K_SYNTH_MAC_INF);
                
                if (!tnfrDynamics && !tnfDepletion)
                {    
                    double il10 = log(((grid.il10(_pos) * MW_IL10 * 1e6)/(NAV * VOL))); // converting il10 concentration to log(ng/mL) for use in dose dependence
                    //std::cout << "LOG OF IL10: " << il10 << std::endl;
                    double tnfMOD = (1.0/(1.0 + exp((il10 + _PARAM(PARAM_GR_LINK_LOG_ALPHA))/_PARAM(PARAM_GR_LINK_LOG_BETA)))); // calculate the fraction of inhibition
                        grid.incTNF(_pos, (tnfMOD * _PARAM(PARAM_MAC_SEC_RATE_TNF) * mdt));
                    //cout << "Debug: IL10 inhibition from Mac::MAC_INFECTED" << std::endl;
                }
                if (!il10rDynamics && !il10Depletion) {
                    grid.incil10(_pos, (_PARAM(PARAM_MAC_SEC_RATE_IL10) * mdt));
                    //cout << "Debug: Secrete from Mac::MAC_INFECTED" << std::endl;
                }
                
                ++grid.nSecretions(_pos);
            }
            
            else
            {
                grid.incCCL2(_pos, (0.5 * _PARAM(PARAM_MAC_SEC_RATE_CCL2) * mdt));
                grid.incCCL5(_pos, (0.5 * _PARAM(PARAM_MAC_SEC_RATE_CCL5) * mdt));
                grid.incCXCL9(_pos,  (0.5 * _PARAM(PARAM_MAC_SEC_RATE_CXCL9) * mdt));
                _kSynth = 0.5 * _PARAM(PARAM_GR_K_SYNTH_MAC);
                _kmRNA = 0.5 * _PARAM(PARAM_GR_K_RNA_MAC);
                _kISynth = _PARAM(PARAM_GR_I_K_SYNTH_MAC_INF);
                
                if (!tnfrDynamics && !tnfDepletion)
                {    
                    double il10 = log(((grid.il10(_pos) * MW_IL10 * 1e6)/(NAV * VOL))); // converting il10 concentration to log(ng/mL) for use in dose dependence
                    //std::cout << "LOG OF IL10: " << il10 << std::endl;
                    double tnfMOD = (1.0/(1.0 + exp((il10 + _PARAM(PARAM_GR_LINK_LOG_ALPHA))/_PARAM(PARAM_GR_LINK_LOG_BETA)))); // calculate the fraction of inhibition
		    grid.incTNF(_pos,(tnfMOD * _PARAM(PARAM_MAC_SEC_RATE_TNF) * mdt));
                    //cout << "Debug: IL10 inhibition from Mac::MAC_INFECTED" << std::endl;
                    //std::cout << "tnfMOD: " << tnfMOD << std::endl;
                }
                if (!il10rDynamics && !il10Depletion) {
                    grid.incil10(_pos, (_PARAM(PARAM_MAC_SEC_RATE_IL10) * mdt));
                    //cout << "Debug: Secrete from Mac::MAC_INFECTED" << std::endl;
                }
                
                ++grid.nSecretions(_pos);
            }
        }
        else if (_state == Mac::MAC_CINFECTED) {
            
            if (_NFkB) {
                grid.incCCL2(_pos, (_PARAM(PARAM_MAC_SEC_RATE_CCL2) * mdt));
                grid.incCCL5(_pos, (_PARAM(PARAM_MAC_SEC_RATE_CCL5) * mdt));
                grid.incCXCL9(_pos,  (_PARAM(PARAM_MAC_SEC_RATE_CXCL9) * mdt));
				_kSynth = _PARAM(PARAM_GR_K_SYNTH_MAC);
                _kmRNA = _PARAM(PARAM_GR_K_RNA_MAC);
                _kISynth = 2.0 * _PARAM(PARAM_GR_I_K_SYNTH_MAC_INF);
                
                if (!tnfrDynamics && !tnfDepletion)
                {    
                    double il10 = log(((grid.il10(_pos) * MW_IL10 * 1e6)/(NAV * VOL))); // converting il10 concentration to log(ng/mL) for use in dose dependence
                    double tnfMOD = (1.0/(1.0 + exp((il10 + _PARAM(PARAM_GR_LINK_LOG_ALPHA))/_PARAM(PARAM_GR_LINK_LOG_BETA)))); // calculate the fraction of inhibition
                        grid.incTNF(_pos, (tnfMOD * _PARAM(PARAM_MAC_SEC_RATE_TNF) * mdt));
                }
                if (!il10rDynamics && !il10Depletion) {
                    grid.incil10(_pos, (2.0 * _PARAM(PARAM_MAC_SEC_RATE_IL10) * mdt));
                }
                
                ++grid.nSecretions(_pos);
            }
        }
        else if (_state == Mac::MAC_ACTIVE) {
            
            if (_NFkB) {
                grid.incCCL2(_pos, (_PARAM(PARAM_MAC_SEC_RATE_CCL2) * mdt));
                grid.incCCL5(_pos, (_PARAM(PARAM_MAC_SEC_RATE_CCL5) * mdt));
                grid.incCXCL9(_pos,  (_PARAM(PARAM_MAC_SEC_RATE_CXCL9) * mdt));
                _kSynth = _PARAM(PARAM_GR_K_SYNTH_MAC);
                _kmRNA = _PARAM(PARAM_GR_K_RNA_MAC);
                _kISynth = _PARAM(PARAM_GR_I_K_SYNTH_MAC_ACT);
                
                if (!tnfrDynamics && !tnfDepletion)
                {    
                    double il10 = log(((grid.il10(_pos) * MW_IL10 * 1e6)/(NAV * VOL))); // converting il10 concentration to log(ng/mL) for use in dose dependence
                    double tnfMOD = (1.0/(1.0 + exp((il10 + _PARAM(PARAM_GR_LINK_LOG_ALPHA))/_PARAM(PARAM_GR_LINK_LOG_BETA)))); // calculate the fraction of inhibition
		    grid.incTNF(_pos, (tnfMOD * _PARAM(PARAM_MAC_SEC_RATE_TNF) * mdt));
                }
                if (!il10rDynamics && !il10Depletion) {
                    //grid.il10(_pos) += (0.0);
                }
                
                ++grid.nSecretions(_pos);
            }
        }
	}
}


void Mac::updateState()
{
	_state = _nextState;
}

void Mac::kill()
{
	_nextState = _state = Mac::MAC_DEAD;
}

void Mac::apoptosis(GrGrid& grid)
{
	_nextState = Mac::MAC_DEAD;

	if (_state == Mac::MAC_INFECTED || _state == Mac::MAC_CINFECTED)
	{
		// disperse bacteria to the Moore neighborhood
		// (half of the intracellular bacteria dies, the other half is dispersed)
		disperseMtb(grid, 0.5);
		_intMtb = 0;
	}
}

void Mac::computeNextState(const int time, GrGrid& grid, GrStat& stats, bool tnfrDynamics, bool nfkbDynamics, bool, bool)
{
	double tnfBoundFraction = grid.TNF(_pos) / (grid.TNF(_pos) + _PARAM(PARAM_GR_KD1) * 48.16e11);
	double nfkb_adjusted_k_apoptosis = _PARAM(PARAM_GR_K_APOPTOSIS_NFKB_MOLECULAR) * (_PARAM(PARAM_GR_K_IAP)/(_PARAM(PARAM_GR_K_IAP) + _normalizedIAP));
	
	
	// check if it is time to die
	if (timeToDie(time))
	{
		// release intracellular bacteria
		grid.extMTB(_pos) += (_intMtb);
		_intMtb = 0;

		// in case active, death contributes to caseation
		if (_state == Mac::MAC_ACTIVE)
		{
			grid.incKillings(_pos);
		}

		_nextState = Mac::MAC_DEAD;
	}
	else if (!nfkbDynamics && tnfrDynamics && _intBoundTNFR1 > _PARAM(PARAM_GR_THRESHOLD_APOPTOSIS_TNF_MOLECULAR) &&	
			 g_Rand.getReal() < 1 - pow(2.7183, -_PARAM(PARAM_GR_K_APOPTOSIS_MOLECULAR) * (_intBoundTNFR1 - _PARAM(PARAM_GR_THRESHOLD_APOPTOSIS_TNF_MOLECULAR))))
	{
		// TNF induced apoptosis
		stats.incMacApoptosisTNF(_state);
		apoptosis(grid);
	}
	else if (nfkbDynamics && _intBoundTNFR1 > _PARAM(PARAM_GR_THRESHOLD_APOPTOSIS_TNF_MOLECULAR) &&	
			 g_Rand.getReal() < 1 - pow(2.7183, -nfkb_adjusted_k_apoptosis * (_intBoundTNFR1 - _PARAM(PARAM_GR_THRESHOLD_APOPTOSIS_TNF_MOLECULAR))))
	{
		// TNF induced apoptosis
		stats.incMacApoptosisTNF(_state);
		apoptosis(grid);
	}
	else if (!tnfrDynamics && tnfBoundFraction > _PARAM(PARAM_GR_THRESHOLD_APOPTOSIS_TNF) &&
			 g_Rand.getReal() < 1 - pow(2.7183, -_PARAM(PARAM_GR_K_APOPTOSIS) * (tnfBoundFraction - _PARAM(PARAM_GR_THRESHOLD_APOPTOSIS_TNF))))
	{
		// TNF induced apoptosis
		stats.incMacApoptosisTNF(_state);
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
			bool tnfInducedNFkB;
			if (tnfrDynamics)
				tnfInducedNFkB = _surfBoundTNFR1 > _PARAM(PARAM_MAC_THRESHOLD_NFKB_TNF_MOLECULAR) && 
				g_Rand.getReal() < 1 - pow(2.7183, -_PARAM(PARAM_MAC_K_NFKB_MOLECULAR) * (_surfBoundTNFR1 - _PARAM(PARAM_MAC_THRESHOLD_NFKB_TNF_MOLECULAR)));
			else
				tnfInducedNFkB = tnfBoundFraction > _PARAM(PARAM_MAC_THRESHOLD_NFKB_TNF) && 
					g_Rand.getReal() < 1 - pow(2.7183, -_PARAM(PARAM_MAC_K_NFKB) * (tnfBoundFraction - _PARAM(PARAM_MAC_THRESHOLD_NFKB_TNF)));
			
			_NFkB = _state == Mac::MAC_CINFECTED || _state == Mac::MAC_ACTIVE || tnfInducedNFkB ||
				getExtMtbInMoore(grid) > _PARAM(PARAM_MAC_THRESHOLD_NFKB_EXTMTB);
				//grid.extMTB(_pos) > _PARAM(PARAM_MAC_THRESHOLD_NFKB_EXTMTB);
			
			if (_stat1 && _surfBoundIL10R < _PARAM(PARAM_MAC_THRESHOLD_ICOS))
			{
				_ICOS = _stat1;
			}

			switch (_state)
			{
			case Mac::MAC_DEAD:
				// if dead, stay dead
				_nextState = Mac::MAC_DEAD;
				break;
			case Mac::MAC_RESTING:
				handleResting(time, grid, stats, nfkbDynamics);
				break;
			case Mac::MAC_INFECTED:
				handleInfected(time, grid, stats, nfkbDynamics);
				break;
			case Mac::MAC_CINFECTED:
				handleChronicallyInfected(time, grid, stats);
				break;
			case Mac::MAC_ACTIVE:
				handleActivated(time, grid, stats);
				break;
      default:
        throw std::runtime_error("Invalid Macrophage State");
			}
		}
	}
}

void Mac::solveDegradation(GrGrid& grid, double dt, bool tnfrDynamics, bool il10rDynamics)
{
	// DBG
//    cout << "Debug: Running T cell degradation"
//    	<<  " _PARAM(PARAM_GR_MEAN_TNFR1_MAC): " << _PARAM(PARAM_GR_MEAN_TNFR1_MAC)
//    	<<  " _PARAM(PARAM_GR_I_IL10R_MAC): " << _PARAM(PARAM_GR_I_IL10R_MAC)
//    	<< std::endl;
    // DBG

	Agent::solveDegradation(grid, dt, tnfrDynamics, il10rDynamics, _PARAM(PARAM_GR_MEAN_TNFR1_MAC), _PARAM(PARAM_GR_I_IL10R_MAC));
}

void Mac::handleResting(const int time, GrGrid& grid, GrStat& stats, bool nfkbDynamics)
{
	_stat1 |= g_Rand.getReal() < getCountTgam(Tgam::TGAM_ACTIVE, grid) * _PARAM(PARAM_MAC_PROB_STAT1_TGAM);

	/* You can get infected from only 1 bacteria;
	 * there should always be a probability associated with getting infected.
	 *
	 * Suggestion:
	 * - always engulf as many bacteria as specified by _EXT_MTB_ENGULF_MAC_REST
	 *   (another suggestion: we can even fix this parameter to 1 and thus get rid of it)
	 * - with a certain probability kill those bacteria, otherwise become infected */
	if (grid.extMTB(_pos) <= _PARAM(PARAM_MAC_NR_UPTAKE_RI_EXTMTB))
	{
		// kill extracellular bacteria, only if there are not too many
		grid.extMTB(_pos) = (0);
		_nextState = Mac::MAC_RESTING;
	}
	else if (g_Rand.getReal() < _PARAM(PARAM_MAC_PROB_KILL_R_EXTMTB))
	{
		// there are too many extracellular bacteria, but we can kill some of those
		// with probability PARAM_MAC_PROB_KILL_R_EXTMTB
		grid.extMTB(_pos) += (-1 * std::min(_PARAM(PARAM_MAC_NR_UPTAKE_RI_EXTMTB), grid.extMTB(_pos)));
		_nextState = Mac::MAC_RESTING;
	}
	else
	{
		// too many bacteria and no killing => macrophage may become infected
		_intMtb = _PARAM(PARAM_MAC_NR_UPTAKE_RI_EXTMTB);
		grid.extMTB(_pos) = (std::max<double>(grid.extMTB(_pos) - _PARAM(PARAM_MAC_NR_UPTAKE_RI_EXTMTB), 0.0));

		_nextState = Mac::MAC_INFECTED;
	}

	// macrophage may become activated
	if (nfkbDynamics) // macrophage activation dynamics follow NFkB dynamics
	{
		if (_stat1)
		{
			if (getExtMtbInMoore(grid) > _PARAM(PARAM_MAC_THRESHOLD_NFKB_EXTMTB))
			{
				_intMtb = 0;
				_activationTime = time;
				_deathTime = time + _PARAM(PARAM_MAC_A_AGE);
				_nextState = Mac::MAC_ACTIVE;
			}
			else if (_normalizedACT > _PARAM(PARAM_GR_ACT_THRESHOLD) &&
					 g_Rand.getReal() < 1 - pow(2.7183, -_PARAM(PARAM_GR_ACT_K) * (_normalizedACT - _PARAM(PARAM_GR_ACT_THRESHOLD))))
			{
				_intMtb = 0;
				_activationTime = time;
				_deathTime = time + _PARAM(PARAM_MAC_A_AGE);
				_nextState = Mac::MAC_ACTIVE;
				stats.incRestingMacActivationTNF();
			}
		}
	}
	else if (_stat1 && _NFkB)
	{
		_intMtb = 0;
		_activationTime = time;
		_deathTime = time + _PARAM(PARAM_MAC_A_AGE);
		_nextState = Mac::MAC_ACTIVE;
		stats.incRestingMacActivationTNF();
	}
}

void Mac::handleInfected(const int time, GrGrid& grid, GrStat& stats, bool nfkbDynamics)
{
	// intracellular bacteria reproduce
	_intMtb *= getIntMtbGrowthRate(time);

	// uptake extracellular bacteria
	// The probability to uptake decreases linearly as the number of intra-cellular bacteria increases to the threshold
	// for becoming chronically infected. At that threshold the probability becomes 0.
	// The probability of uptake is the compliment of the probability of killing.
	double baseProbExtMtbUptake = (1.0 - _PARAM(PARAM_MAC_PROB_KILL_R_EXTMTB))/2.0;
	double probExtMtbUptake = (baseProbExtMtbUptake *  (1.0 - (_intMtb / _PARAM(PARAM_MAC_THRESHOLD_BECOME_CI_INTMTB)) ));
	if (grid.extMTB(_pos) > 0 && g_Rand.getReal() < probExtMtbUptake)
	{
		double dExtMtb = std::min(grid.extMTB(_pos), _PARAM(PARAM_MAC_NR_UPTAKE_RI_EXTMTB));
		grid.extMTB(_pos) += (-dExtMtb);
		_intMtb += dExtMtb;
	}

	if (_intMtb >= _PARAM(PARAM_MAC_THRESHOLD_BECOME_CI_INTMTB))
	{
		_nextState = Mac::MAC_CINFECTED;
	}
	else
	{
		_stat1 |= g_Rand.getReal() <
			getCountTgam(Tgam::TGAM_ACTIVE, grid) * _PARAM(PARAM_MAC_PROB_STAT1_TGAM);

		// macrophage may become activated
		if (nfkbDynamics) // macrophage activation dynamics follow NFkB dynamics
		{
			if (_stat1)
			{
				if (getExtMtbInMoore(grid) > _PARAM(PARAM_MAC_THRESHOLD_NFKB_EXTMTB))
				{
					_intMtb = 0;
					_activationTime = time;
					_deathTime = time + _PARAM(PARAM_MAC_A_AGE);
					_nextState = Mac::MAC_ACTIVE;
				}
				else if (_normalizedACT > _PARAM(PARAM_GR_ACT_THRESHOLD) &&
						 g_Rand.getReal() < 1 - pow(2.7183, -_PARAM(PARAM_GR_ACT_K) * (_normalizedACT - _PARAM(PARAM_GR_ACT_THRESHOLD))))
				{
					_intMtb = 0;
					_activationTime = time;
					_deathTime = time + _PARAM(PARAM_MAC_A_AGE);
					_nextState = Mac::MAC_ACTIVE;
					stats.incInfMacActivationTNF();
				}
				else 
				{
					_nextState = Mac::MAC_INFECTED;
				}
			}
		}
		else 
		{
			if (_stat1 && _NFkB)
			{
				_intMtb = 0;
				_activationTime = time;
				_deathTime = time + _PARAM(PARAM_MAC_A_AGE);
				_nextState = Mac::MAC_ACTIVE;
				stats.incInfMacActivationTNF();
			}
			else
			{
				_nextState = Mac::MAC_INFECTED;
			}
		}
	}
}

void Mac::handleChronicallyInfected(const int time, GrGrid& grid, GrStat&)
{
	// intracellular bacteria reproduce
	_intMtb *= getIntMtbGrowthRate(time);

	if (_intMtb >= _PARAM(PARAM_MAC_THRESHOLD_BURST_CI_INTMTB))
	{
		// burst, all intracellular bacteria disperse to the Moore neighborhood
		disperseMtb(grid, 1.0);
		_intMtb = 0;

		// increment number of killings
		grid.incKillings(_pos);

		_nextState = Mac::MAC_DEAD;
	}
	else
	{
		_nextState = Mac::MAC_CINFECTED;
	}
}

void Mac::handleActivated(const int, GrGrid& grid, GrStat&)
{
	// kill extracellular bacteria in the compartment the macrophage resides
	grid.extMTB(_pos) += (-1 * std::min(grid.extMTB(_pos), _PARAM(PARAM_MAC_NR_UPTAKE_A_EXTMTB)));

	_nextState = Mac::MAC_ACTIVE;
}

void Mac::deactivate(const int time)
{
	switch (_state)
	{
	case Mac::MAC_RESTING:
	case Mac::MAC_INFECTED:
	case Mac::MAC_ACTIVE:
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
	case Mac::MAC_DEAD:
		std::cout << "dead, ";
		break;
	case Mac::MAC_RESTING:
		std::cout << "resting, ";
		break;
	case Mac::MAC_INFECTED:
		std::cout << "infected, ";
		break;
	case Mac::MAC_CINFECTED:
		std::cout << "chronically infected, ";
		break;
	case Mac::MAC_ACTIVE:
		std::cout << "active, ";
		break;
  default: break;
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

	Serialization::writeHeader(out, Mac::_ClassName);

	Agent::serialize(out);

    out << _macodeSize << std::endl;
    
	int intVal = (int) _state;
	out << intVal << std::endl;

	intVal = (int) _nextState;
	out << intVal << std::endl;

	out << _intMtb << std::endl;
	out << _NFkB << std::endl;
	out << _stat1 << std::endl;
	out << _ICOS << std::endl;
	out << _activationTime << std::endl;
	out << _deactivationTime << std::endl;

	Serialization::writeFooter(out, Mac::_ClassName);
}

void Mac::deserialize(std::istream& in)
{
	assert(in.good());

	if (!Serialization::readHeader(in, Mac::_ClassName))
	{
		exit(1);
	}

	int intVal;

	Agent::deserialize(in);

    in >> _macodeSize;
    
	in >> intVal;
	_state = (Mac::State) intVal;

	in >> intVal;
	_nextState = (Mac::State) intVal;

	in >> _intMtb;
	in >> _NFkB;
	in >> _stat1;
	in >> _ICOS;
	in >> _activationTime;
	in >> _deactivationTime;
	
	if (!Serialization::readFooter(in, Mac::_ClassName))
	{
		exit(1);
	}
}

int Mac::getCountTgam(Tgam::State state, const GrGrid& grid) const
{
	int count = 0;

	// count the number of Tgam cells in the Moore neighborhood
	// that are in the specified state
	for (int i = -1; i <= 1; i++)
		for (int j = -1; j <= 1; j++)
    {
      Pos p(grid.mod_row(_pos.x + i), grid.mod_col(_pos.y + j));
			count += grid.hasAgentType(TGAM, state, p);
    }

	return count;
}

double Mac::getExtMtbInMoore(const GrGrid& grid) const
{
	double mtb = 0;

	// count the number of extMtb cells in the Moore neighborhood
	// that are in the specified state
	for (int i = -1; i <= 1; i++)
		for (int j = -1; j <= 1; j++) {
      Pos p(grid.mod_row(_pos.x + i), grid.mod_col(_pos.y + j));
      mtb += grid.extMTB(p);
    }

	return mtb;
}

// Disperse the bacteria in a discretized manner. If the macrophage has 1 or more bacteria,
// then disperse to so that each compartment dispersed to receives at least 1 bacterium.
// This means disperse to fewer than all the compartments in the Mooore neighborhood if
// necessary, with the destination compartments chosen randomly. If the number of
// compartments to dispsers to is less than all of the
// For example, if 6.3 bacteria are to be dispersed then disperse to 1.05 bacteria to 6
// randomly chosen neighbors. If 10.4 bacteria are to be dispersed then disperse to 1.16
// bacteria to all compartments in the Moore neighborhood. If 0.8 If 10.4 bacteria are to
// be dispersed then disperse them to 1 compartment.
//
void Mac::disperseMtb(GrGrid& grid, double fraction)
{
    assert(0 < fraction && fraction <= 1);

    std::vector<int> permutation;
    for (int i = 0; i < MOORE_COUNT; i++)
    {
    	permutation.push_back(i);
    }

    double mtb = fraction * _intMtb;

    // Disperse to at least one compartment.
    double minCompartments = std::max(floor(mtb), 1.0);
    double dMtb = mtb / std::min(minCompartments, static_cast<double>(permutation.size()));

    size_t n = std::min(static_cast<size_t>(minCompartments), permutation.size());
    if (n < permutation.size())
        std::random_shuffle(permutation.begin(), permutation.end(), g_Rand);

    for (size_t i = 0; i < n; i++)
    {
    	Pos pos = compartmentOrdinalToCoordinates(permutation[i], grid.getRange());

    	grid.extMTB(pos) += (dMtb);
    }

    _intMtb -= mtb;
}

