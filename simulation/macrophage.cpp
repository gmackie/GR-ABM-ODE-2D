/*
 * macrophage.cpp
 *
 *  Created on: 05-nov-2009
 *      Author: M. El-Kebir
 */

#include "macrophage.h"
#include "params.h"
#include "grgrid.h"
#include "serialization.h"

const std::string Mac::_ClassName = "Mac";

// Needed for deserializing the model state.
// Avoids the calls to the random number generator in the normal constructor, allowing the random number generator
// to remain in synch after deserialization.
Mac::Mac()
	:Agent()
	, _state(MAC_DEAD)
	, _nextState(MAC_DEAD)
	, _intMtb(-1.0)
	, _NFkB(0)
	, _stat1(0)
	, _ICOS(0)
	, _activationTime(-1)
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
    
	// NF-kB signaling pathway components
	, _IKKKa(-1.0)
	, _IKKn(-1.0)
	, _IKKa(-1.0)
	, _IKKi(-1.0)
	, _IkBp(-1.0)
	, _NFkB_IkBp(-1.0)
	, _NFkBc(-1.0)
	, _NFkBn(-1.0)
	, _A20(-1.0)
	, _A20t(-1.0)
	, _IkB(-1.0)
	, _IkBn(-1.0)
	, _IkBt(-1.0)
	, _NFkB_IkB(-1.0)
	, _NFkB_IkBn(-1.0)
	, _GA20(-1.0)
	, _GIkB(-1.0)
	, _GR(-1.0)
	, _c1rrChemTNF(-1.0)
	, _c1rChem(-1.0)
	, _c1rTNF(-1.0)
	, _chemt(-1.0)
	, _chem(-1.0)
	, _TNFt(-1.0)
	, _TNF(-1.0)
	, _ACTt(-1.0)
	, _ACT(-1.0)
	, _normalizedACT(-1.0)
	, _IAPt(-1.0)
	, _IAP(-1.0)
	, _normalizedIAP(-1.0)

{
}

Mac::Mac(int birthtime, int row, int col, MacState state, double intMtb, bool NFkB, bool stat1)
	: Agent(birthtime, birthtime + _PARAM(PARAM_MAC_AGE), row, col)
	, _state(state)
	, _nextState(state)
	, _intMtb(intMtb)
	, _NFkB(NFkB)
	, _stat1(stat1)
	, _ICOS(0)
	, _activationTime(-1)
	, _deactivationTime(-1)
	, _mTNF(0.0)
	, _surfTNFR1(g_Rand.getReal(_PARAM(PARAM_GR_MEAN_TNFR1_MAC)-_PARAM(PARAM_GR_STD_TNFR1_MAC),_PARAM(PARAM_GR_MEAN_TNFR1_MAC)+_PARAM(PARAM_GR_STD_TNFR1_MAC)))
	, _surfTNFR2(g_Rand.getReal(_PARAM(PARAM_GR_MEAN_TNFR2_MAC)-_PARAM(PARAM_GR_STD_TNFR2_MAC),_PARAM(PARAM_GR_MEAN_TNFR2_MAC)+_PARAM(PARAM_GR_STD_TNFR2_MAC)))
//	, _surfTNFR1(g_Rand.getLogNormal(_PARAM(PARAM_GR_MEAN_TNFR1_MAC),_PARAM(PARAM_GR_STD_TNFR1_MAC)))
//	, _surfTNFR2(g_Rand.getLogNormal(_PARAM(PARAM_GR_MEAN_TNFR2_MAC),_PARAM(PARAM_GR_STD_TNFR2_MAC)))
	, _surfBoundTNFR1(0.0)
	, _surfBoundTNFR2(0.0)
	, _intBoundTNFR1(0.0)
	, _intBoundTNFR2(0.0)
    , _mTNFRNA(0.0)
	, _vTNFR1(_surfTNFR1 * _PARAM(PARAM_GR_K_T1))
	, _vTNFR2(_surfTNFR2 * _PARAM(PARAM_GR_K_T2))
	, _kSynth(0.0)
	, _kTACE(_PARAM(PARAM_GR_K_TACE_MAC))
    , _kmRNA(0.0)

    // IL10 components

    , _surfIL10R(g_Rand.getReal(_PARAM(PARAM_GR_I_IL10R_MAC) - _PARAM(PARAM_GR_STD_IL10R_MAC), _PARAM(PARAM_GR_I_IL10R_MAC) + _PARAM(PARAM_GR_STD_IL10R_MAC)))
    , _vIL10R(_surfIL10R * _PARAM(PARAM_GR_I_K_T))
    , _surfBoundIL10R(0.0)
    , _kISynth(0.0)

	// NF-kB signaling pathway components
	, _IKKKa(0.0)
	, _IKKn(_PARAM(PARAM_GR_KNN))
	, _IKKa(0.0)
	, _IKKi(0.0)
	, _IkBp(0.0)
	, _NFkB_IkBp(0.0)
	, _NFkBc(0.0)
	, _NFkBn(0.0)
	, _A20(0.0)
	, _A20t(0.0)
	, _IkB(0.0)
	, _IkBn(0.0)
	, _IkBt(0.0)
	, _NFkB_IkB((_PARAM(PARAM_GR_MEAN_NFKB) > 0 ) ? g_Rand.getLogNormal(_PARAM(PARAM_GR_MEAN_NFKB),0.64872*_PARAM(PARAM_GR_MEAN_NFKB)) : 0)
	, _NFkB_IkBn(0.0)
	, _GA20(0.0)
	, _GIkB(0.0)
	, _GR(0.0)
	, _c1rrChemTNF(0.0)
	, _c1rChem(0.0)
	, _c1rTNF(0.0)
	, _chemt(0.0)
	, _chem(0.0)
	, _TNFt(0.0)
	, _TNF(0.0)
	, _ACTt(0.0)
	, _ACT(0.0)
	, _normalizedACT(0.0)
	, _IAPt(0.0)
	, _IAP(0.0)
	, _normalizedIAP(0.0)

{
}

Mac::~Mac()
{
}

void Mac::move(GrGrid& grid)
{
	Pos pos = Agent::moveAgent(grid, true, true, false, true, _PARAM(PARAM_MAC_MOVEMENT_BONUSFACTOR));

	int newRow = pos.first;
	int newCol = pos.second;

	GridCell& cell = grid(_row, _col);
	GridCell& newCell = grid(newRow, newCol);
	if (newCell.addAgent(this))
	{
		_row = newRow;
		_col = newCol;
		assert_res(cell.removeAgent(this));
	}
}

void Mac::secrete(GrGrid& grid, bool tnfrDynamics, bool nfkbDynamics, bool tnfDepletion, bool il10rDynamics, bool il10Depletion)
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
	GridCell& cell = grid(_row, _col);
	
	if (nfkbDynamics) // TNF and chemokines are secreted as a function of NFkB dynamics
	{
		if (_state == MAC_RESTING)
		{
			_c1rChem = _PARAM(PARAM_GR_epsilon2) * _PARAM(PARAM_GR_c1r);
			_c1rTNF = _PARAM(PARAM_GR_epsilon2) * _PARAM(PARAM_GR_c1r);
			_c1rrChemTNF = 0;
            _kISynth = 0.0;
		}
		else if (_state == MAC_INFECTED)
		{
			_c1rChem = _PARAM(PARAM_GR_c1r);
			_c1rTNF = _PARAM(PARAM_GR_c1r);
			_c1rrChemTNF = _PARAM(PARAM_GR_epsilon1) * _PARAM(PARAM_GR_c1r);
			cell.incNrSecretions();
            _kISynth = _PARAM(PARAM_GR_I_K_SYNTH_MAC_INF);
            
            if (!il10rDynamics && !il10Depletion) {
                cell.incIL10(_PARAM(PARAM_MAC_SEC_RATE_IL10));
            }
            
		}
		else if (_state == MAC_ACTIVE)
		{
			_c1rChem = _PARAM(PARAM_GR_c1r);
			_c1rTNF = _PARAM(PARAM_GR_c1r);
			_c1rrChemTNF = _PARAM(PARAM_GR_c1r);
			cell.incNrSecretions();
            _kISynth = _PARAM(PARAM_GR_I_K_SYNTH_MAC_ACT);
            
            if (!il10rDynamics && !il10Depletion) {
                cell.incIL10(0.5 * _PARAM(PARAM_MAC_SEC_RATE_IL10));
            }
            
		}
        else if (_state == MAC_CINFECTED)
        {
            _c1rChem = _PARAM(PARAM_GR_c1r);
			_c1rTNF = _PARAM(PARAM_GR_c1r);
			_c1rrChemTNF = _PARAM(PARAM_GR_c1r);
			cell.incNrSecretions();
            _kISynth = 2.0 * _PARAM(PARAM_GR_I_K_SYNTH_MAC_INF);
            
            if (!il10rDynamics && !il10Depletion) {
                cell.incIL10(2.0 * _PARAM(PARAM_MAC_SEC_RATE_IL10));
            }
        }
	}
    
	else
	// TNF, IL10, and chemokines are secreted independent of NFkB dynamics
    // IL10 dynamics are independent of the state of _NFkB
	{
        double Nav = 6.02e23; // Avogadro Number
        double vol = 8.0e-12; // volume of a cell in liter
        double MW_IL10 = 18600; // molecular weight of IL10 in g/mol
        
		if (_state == MAC_RESTING) {
            
            if (_NFkB) {
                cell.incCCL2(0.5 * _PARAM(PARAM_MAC_SEC_RATE_CCL2));
				cell.incCCL5(0.5 * _PARAM(PARAM_MAC_SEC_RATE_CCL5));
				cell.incCXCL9(0.5 * _PARAM(PARAM_MAC_SEC_RATE_CXCL9));
				_kSynth = 0.5 * _PARAM(PARAM_GR_K_SYNTH_MAC);
                _kmRNA = 0.5 * _PARAM(PARAM_GR_K_RNA_MAC);
                _kISynth = 0.0;
                
                if (!tnfrDynamics && !tnfDepletion)
                {
                    
                    double il10 = log(((cell.getIL10() * MW_IL10 * 1e6)/(Nav * vol))); // converting il10 concentration to log(ng/mL) for use in dose dependence
                    //std::cout << "LOG OF IL10: " << il10 << std::endl;
                    double tnfMOD = (1.0/(1.0 + exp((il10 + _PARAM(PARAM_GR_LINK_LOG_ALPHA))/_PARAM(PARAM_GR_LINK_LOG_BETA)))); // calculate the fraction of inhibition
                    //std::cout << "tnfMOD: " << tnfMOD << std::endl;
					cell.incTNF(tnfMOD * 0.5 * _PARAM(PARAM_MAC_SEC_RATE_TNF));
                }
                
                cell.incNrSecretions();
            }
            else
            {
                _kSynth = 0.0;
                _kmRNA = 0.0;
                _kISynth = 0.0;
            }
        }
        else if (_state == MAC_INFECTED) {
            
            if (_NFkB) {
                cell.incCCL2(_PARAM(PARAM_MAC_SEC_RATE_CCL2));
				cell.incCCL5(_PARAM(PARAM_MAC_SEC_RATE_CCL5));
				cell.incCXCL9(_PARAM(PARAM_MAC_SEC_RATE_CXCL9));
				_kSynth = _PARAM(PARAM_GR_K_SYNTH_MAC);
                _kmRNA = _PARAM(PARAM_GR_K_RNA_MAC);
                _kISynth = _PARAM(PARAM_GR_I_K_SYNTH_MAC_INF);
                
                if (!tnfrDynamics && !tnfDepletion)
                {    
                    double il10 = log(((cell.getIL10() * MW_IL10 * 1e6)/(Nav * vol))); // converting il10 concentration to log(ng/mL) for use in dose dependence
                    //std::cout << "LOG OF IL10: " << il10 << std::endl;
                    double tnfMOD = (1.0/(1.0 + exp((il10 + _PARAM(PARAM_GR_LINK_LOG_ALPHA))/_PARAM(PARAM_GR_LINK_LOG_BETA)))); // calculate the fraction of inhibition
					cell.incTNF(tnfMOD * _PARAM(PARAM_MAC_SEC_RATE_TNF));
                    //cout << "Debug: IL10 inhibition from MAC_INFECTED" << std::endl;
                }
                if (!il10rDynamics && !il10Depletion) {
                    cell.incIL10(_PARAM(PARAM_MAC_SEC_RATE_IL10));
                    //cout << "Debug: Secrete from MAC_INFECTED" << std::endl;
                }
                
                cell.incNrSecretions();
            }
            
            else
            {
                cell.incCCL2(0.5 * _PARAM(PARAM_MAC_SEC_RATE_CCL2));
                cell.incCCL5(0.5 * _PARAM(PARAM_MAC_SEC_RATE_CCL5));
                cell.incCXCL9(0.5 * _PARAM(PARAM_MAC_SEC_RATE_CXCL9));
                _kSynth = 0.5 * _PARAM(PARAM_GR_K_SYNTH_MAC);
                _kmRNA = 0.5 * _PARAM(PARAM_GR_K_RNA_MAC);
                _kISynth = _PARAM(PARAM_GR_I_K_SYNTH_MAC_INF);
                
                if (!tnfrDynamics && !tnfDepletion)
                {    
                    double il10 = log(((cell.getIL10() * MW_IL10 * 1e6)/(Nav * vol))); // converting il10 concentration to log(ng/mL) for use in dose dependence
                    //std::cout << "LOG OF IL10: " << il10 << std::endl;
                    double tnfMOD = (1.0/(1.0 + exp((il10 + _PARAM(PARAM_GR_LINK_LOG_ALPHA))/_PARAM(PARAM_GR_LINK_LOG_BETA)))); // calculate the fraction of inhibition
					cell.incTNF(tnfMOD * _PARAM(PARAM_MAC_SEC_RATE_TNF));
                    //cout << "Debug: IL10 inhibition from MAC_INFECTED" << std::endl;
                    //std::cout << "tnfMOD: " << tnfMOD << std::endl;
                }
                if (!il10rDynamics && !il10Depletion) {
                    cell.incIL10(_PARAM(PARAM_MAC_SEC_RATE_IL10));
                    //cout << "Debug: Secrete from MAC_INFECTED" << std::endl;
                }
                
                cell.incNrSecretions();
            }
        }
        else if (_state == MAC_CINFECTED) {
            
            if (_NFkB) {
                cell.incCCL2(_PARAM(PARAM_MAC_SEC_RATE_CCL2));
				cell.incCCL5(_PARAM(PARAM_MAC_SEC_RATE_CCL5));
				cell.incCXCL9(_PARAM(PARAM_MAC_SEC_RATE_CXCL9));
				_kSynth = _PARAM(PARAM_GR_K_SYNTH_MAC);
                _kmRNA = _PARAM(PARAM_GR_K_RNA_MAC);
                _kISynth = _kISynth = 2.0 * _PARAM(PARAM_GR_I_K_SYNTH_MAC_INF);
                
                if (!tnfrDynamics && !tnfDepletion)
                {    
                    double il10 = log(((cell.getIL10() * MW_IL10 * 1e6)/(Nav * vol))); // converting il10 concentration to log(ng/mL) for use in dose dependence
                    double tnfMOD = (1.0/(1.0 + exp((il10 + _PARAM(PARAM_GR_LINK_LOG_ALPHA))/_PARAM(PARAM_GR_LINK_LOG_BETA)))); // calculate the fraction of inhibition
					cell.incTNF(tnfMOD * _PARAM(PARAM_MAC_SEC_RATE_TNF));
                }
                if (!il10rDynamics && !il10Depletion) {
                    cell.incIL10(2.0 * _PARAM(PARAM_MAC_SEC_RATE_IL10));
                }
                
                cell.incNrSecretions();
            }
        }
        else if (_state == MAC_ACTIVE) {
            
            if (_NFkB) {
                cell.incCCL2(_PARAM(PARAM_MAC_SEC_RATE_CCL2));
				cell.incCCL5(_PARAM(PARAM_MAC_SEC_RATE_CCL5));
				cell.incCXCL9(_PARAM(PARAM_MAC_SEC_RATE_CXCL9));
				_kSynth = _PARAM(PARAM_GR_K_SYNTH_MAC);
                _kmRNA = _PARAM(PARAM_GR_K_RNA_MAC);
                _kISynth = _PARAM(PARAM_GR_I_K_SYNTH_MAC_ACT);
                
                if (!tnfrDynamics && !tnfDepletion)
                {    
                    double il10 = log(((cell.getIL10() * MW_IL10 * 1e6)/(Nav * vol))); // converting il10 concentration to log(ng/mL) for use in dose dependence
                    double tnfMOD = (1.0/(1.0 + exp((il10 + _PARAM(PARAM_GR_LINK_LOG_ALPHA))/_PARAM(PARAM_GR_LINK_LOG_BETA)))); // calculate the fraction of inhibition
					cell.incTNF(tnfMOD * _PARAM(PARAM_MAC_SEC_RATE_TNF));
                }
                if (!il10rDynamics && !il10Depletion) {
                    cell.incIL10(0.0);
                }
                
                cell.incNrSecretions();
            }
        }
        
        
//        if (_NFkB)
//		{
//			// if NFkB is turned on, secrete TNF and chemokines
//			if (_state != MAC_RESTING)
//			{
//				cell.incCCL2(_PARAM(PARAM_MAC_SEC_RATE_CCL2));
//				cell.incCCL5(_PARAM(PARAM_MAC_SEC_RATE_CCL5));
//				cell.incCXCL9(_PARAM(PARAM_MAC_SEC_RATE_CXCL9));
//				_kSynth = _PARAM(PARAM_GR_K_SYNTH_MAC);
//                
//                if (_state == MAC_ACTIVE)
//                {
//                    _kISynth = _PARAM(PARAM_GR_I_K_SYNTH_MAC_ACT);
//                }
//                if (_state == MAC_INFECTED)
//                {
//                    _kISynth = _PARAM(PARAM_GR_I_K_SYNTH_MAC_INF);
//                }
//                if (_state == MAC_CINFECTED)
//                {
//                _kISynth = _kISynth = 2.0 * _PARAM(PARAM_GR_I_K_SYNTH_MAC_INF);
//                }
//                                      
//				if (!tnfrDynamics && !tnfDepletion)
//					cell.incTNF(_PARAM(PARAM_MAC_SEC_RATE_TNF));
//                
//			}
//			else
//			{
//				cell.incCCL2(0.5 * _PARAM(PARAM_MAC_SEC_RATE_CCL2));
//				cell.incCCL5(0.5 * _PARAM(PARAM_MAC_SEC_RATE_CCL5));
//				cell.incCXCL9(0.5 * _PARAM(PARAM_MAC_SEC_RATE_CXCL9));
//				_kSynth = 0.5 * _PARAM(PARAM_GR_K_SYNTH_MAC);
//                _kISynth = 0.0;
//				if (!tnfrDynamics && !tnfDepletion)
//					cell.incTNF(0.5 * _PARAM(PARAM_MAC_SEC_RATE_TNF));
//			}
//		
//			cell.incNrSecretions();
//		}
//		else if (_state == MAC_RESTING)
//        {
//			_kSynth = 0.0;
//            _kISynth = 0.0;
//        }
//		else if (_state == MAC_INFECTED)
//		{
//			cell.incCCL2(0.5 * _PARAM(PARAM_MAC_SEC_RATE_CCL2));
//			cell.incCCL5(0.5 * _PARAM(PARAM_MAC_SEC_RATE_CCL5));
//			cell.incCXCL9(0.5 * _PARAM(PARAM_MAC_SEC_RATE_CXCL9));
//			_kSynth = 0.5 * _PARAM(PARAM_GR_K_SYNTH_MAC);
//            _kISynth = _PARAM(PARAM_GR_I_K_SYNTH_MAC_INF);
//			if (!tnfrDynamics && !tnfDepletion)
//				cell.incTNF(0.5 * _PARAM(PARAM_MAC_SEC_RATE_TNF));
//		
//			cell.incNrSecretions();
//		}
//        else if (_state == MAC_CINFECTED)
//        {
//            _kISynth = _kISynth = 2.0 * _PARAM(PARAM_GR_I_K_SYNTH_MAC_INF);
//        }
//        
//        else if (_state == MAC_ACTIVE)
//        {
//            _kISynth = _PARAM(PARAM_GR_I_K_SYNTH_MAC_ACT);
//        }
        
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

	if (_state == MAC_INFECTED || _state == MAC_CINFECTED)
	{
		// disperse bacteria to the Moore neighborhood
		// (half of the intracellular bacteria dies, the other half is dispersed)
		disperseMtb(grid, 0.5);
		_intMtb = 0;
	}
}

void Mac::computeNextState(const int time, GrGrid& grid, GrStat& stats, bool tnfrDynamics, bool nfkbDynamics, bool)
{
	GridCell& cell = grid(_row, _col);
	double tnfBoundFraction = cell.getTNF() / (cell.getTNF() + _PARAM(PARAM_GR_KD1) * 48.16e11);
	double nfkb_adjusted_k_apoptosis = _PARAM(PARAM_GR_K_APOPTOSIS_NFKB_MOLECULAR) * (_PARAM(PARAM_GR_K_IAP)/(_PARAM(PARAM_GR_K_IAP) + _normalizedIAP));
	
	
	// check if it is time to die
	if (timeToDie(time))
	{
		// release intracellular bacteria
		cell.incExtMtb(_intMtb);
		_intMtb = 0;

		// in case active, death contributes to caseation
		if (_state == MAC_ACTIVE)
		{
			cell.incNrKillings();
		}

		_nextState = MAC_DEAD;
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
			
			_NFkB = _state == MAC_CINFECTED || _state == MAC_ACTIVE || tnfInducedNFkB ||
				getExtMtbInMoore(grid) > _PARAM(PARAM_MAC_THRESHOLD_NFKB_EXTMTB);
				//cell.getExtMtb() > _PARAM(PARAM_MAC_THRESHOLD_NFKB_EXTMTB);
			
			if (_stat1 && _surfBoundIL10R < _PARAM(PARAM_MAC_THRESHOLD_ICOS))
			{
				_ICOS = _stat1;
			}

			switch (_state)
			{
			case MAC_DEAD:
				// if dead, stay dead
				_nextState = MAC_DEAD;
				break;
			case MAC_RESTING:
				handleResting(time, grid, stats, nfkbDynamics);
				break;
			case MAC_INFECTED:
				handleInfected(time, grid, stats, nfkbDynamics);
				break;
			case MAC_CINFECTED:
				handleChronicallyInfected(time, grid, stats);
				break;
			case MAC_ACTIVE:
				handleActivated(time, grid, stats);
				break;
      default:
        throw std::runtime_error("Invalid Macrophage State");
			}
		}
	}
}

void Mac::solveTNF(GrGrid& grid, double dt)
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
    
    //cout << "Debug: Running TNF dynamics" << std::endl;
}

void Mac::solveNFkBandTNF(GrGrid& grid, double dt)
{
	GridCell& cell = grid(_row, _col);
	
	double koff1 = _PARAM(PARAM_GR_K_ON1) * _PARAM(PARAM_GR_KD1);
	double koff2 = _PARAM(PARAM_GR_K_ON2) * _PARAM(PARAM_GR_KD2);
	double kdeg = 4.58e-4;
	double density = 1.25e11; // used for conversion of conc. unit (M -> #/cell) based on cell and microcompartment volumes 
	double Nav = 6.02e23; // Avogadro Number
	double vol = 8.0e-12; // volume of a cell in liter
	double kv = 5; //ratio of cytoplasmic to nuclear volume
	
	double tnf = cell.getTNF() / (Nav * vol);
	double shedtnfr2 = cell.getShedTNFR2() / (Nav * vol);
	double il10 = cell.getIL10() /(Nav * vol);
    
	double dmTNF;
	double dsurfTNFR1;
	double dsurfTNFR2;
	double dsurfBoundTNFR1;
	double dsurfBoundTNFR2;
	double dintBoundTNFR1;
	double dintBoundTNFR2;
	double dsTNF;
	double dshedTNFR2;
	
	double dIKKKa; 
	double dIKKn; 
	double dIKKa; 
	double dIKKi; 
	double dIkBp; 
	double dNFkB_IkBp; 
	double dNFkBc; 
	double dNFkBn; 
	double dA20;
	double dA20t; 
	double dIkB;
	double dIkBn; 
	double dIkBt; 
	double dNFkB_IkB;
	double dNFkB_IkBn;
	double dGA20;
	double dGIkB;
	double dGR; 
	double dchemt; 
	double dchem; 
	double dTNFt;
	double dTNF; 
	double dACTt; 
	double dACT; 
	double dIAPt; 
	double dIAP;
    
    
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
    
    
    // CURRENTLY NOT LINKED WITH IL10 SINCE NFKB DYNAMICS HAVE NOT BEEN LOOKED AT
    
    // tnf dynamics model equations
	dmTNF = (_PARAM(PARAM_GR_e3TNF)*_TNF - _kTACE * _mTNF) * dt;
	dsurfTNFR1 = (_vTNFR1 - _PARAM(PARAM_GR_K_ON1) * tnf * _surfTNFR1 + (koff1+kdeg) * _surfBoundTNFR1 - _PARAM(PARAM_GR_K_T1) * _surfTNFR1 + _PARAM(PARAM_GR_K_REC1) * _intBoundTNFR1) * dt;
	dsurfTNFR2 = (_vTNFR2 - _PARAM(PARAM_GR_K_ON2) * tnf * _surfTNFR2 + (koff2+kdeg) * _surfBoundTNFR2 - _PARAM(PARAM_GR_K_T2) * _surfTNFR2 + _PARAM(PARAM_GR_K_REC2) * _intBoundTNFR2) * dt;
	dsurfBoundTNFR1 = (_PARAM(PARAM_GR_K_ON1) * tnf * _surfTNFR1 - (koff1+kdeg) * _surfBoundTNFR1 - _PARAM(PARAM_GR_K_INT1) * _surfBoundTNFR1) * dt;
	dsurfBoundTNFR2 = (_PARAM(PARAM_GR_K_ON2) * tnf * _surfTNFR2 - (koff2+kdeg) * _surfBoundTNFR2 - _PARAM(PARAM_GR_K_INT2) * _surfBoundTNFR2 - _PARAM(PARAM_GR_K_SHED) * _surfBoundTNFR2) * dt;
	dintBoundTNFR1 = (_PARAM(PARAM_GR_K_INT1) * _surfBoundTNFR1 - _PARAM(PARAM_GR_K_DEG1) * _intBoundTNFR1 - _PARAM(PARAM_GR_K_REC1) * _intBoundTNFR1) * dt;
	dintBoundTNFR2 = (_PARAM(PARAM_GR_K_INT2) * _surfBoundTNFR2 - _PARAM(PARAM_GR_K_DEG2) * _intBoundTNFR2 - _PARAM(PARAM_GR_K_REC2) * _intBoundTNFR2) * dt;
	dsTNF = ((density/Nav) * (_kTACE * _mTNF - _PARAM(PARAM_GR_K_ON1) * tnf * _surfTNFR1 + koff1 * _surfBoundTNFR1 - _PARAM(PARAM_GR_K_ON2) * tnf * _surfTNFR2 + koff2 * _surfBoundTNFR2)) * dt; 
	dshedTNFR2 = ((density/Nav) * _PARAM(PARAM_GR_K_SHED) * _surfBoundTNFR2) * dt;
	
	// NF-kB dynamics model equations
	dIKKKa = (_PARAM(PARAM_GR_ka)*_surfBoundTNFR1*(_PARAM(PARAM_GR_KN)-_IKKKa)*_PARAM(PARAM_GR_kA20)/(_PARAM(PARAM_GR_kA20)+_A20)-_PARAM(PARAM_GR_ki)*_IKKKa) * dt;
	dIKKn = (_PARAM(PARAM_GR_k4)*(_PARAM(PARAM_GR_KNN)-_IKKn-_IKKa-_IKKi)-_PARAM(PARAM_GR_k1)*pow(_IKKKa,2.0)*_IKKn) * dt;
	dIKKa = (_PARAM(PARAM_GR_k1)*pow(_IKKKa,2.0)*_IKKn-_PARAM(PARAM_GR_k3)*_IKKa*(_PARAM(PARAM_GR_k2)+_A20)/_PARAM(PARAM_GR_k2)) * dt;
	dIKKi = (_PARAM(PARAM_GR_k3)*_IKKa*(_PARAM(PARAM_GR_k2)+_A20)/_PARAM(PARAM_GR_k2)-_PARAM(PARAM_GR_k4)*_IKKi) * dt;
	dIkBp = (_PARAM(PARAM_GR_a2)*_IKKa*_IkB-_PARAM(PARAM_GR_tp)*_IkBp) * dt;
	dNFkB_IkBp = (_PARAM(PARAM_GR_a3)*_IKKa*_NFkB_IkB-_PARAM(PARAM_GR_tp)*_NFkB_IkBp) * dt;
	dNFkBc = (_PARAM(PARAM_GR_c6a)*_NFkB_IkB-_PARAM(PARAM_GR_a1)*_NFkBc*_IkB+_PARAM(PARAM_GR_tp)*_NFkB_IkBp-_PARAM(PARAM_GR_i1)*_NFkBc) * dt;
	dNFkBn = (_PARAM(PARAM_GR_i1)*_NFkBc-_PARAM(PARAM_GR_a1)*kv*_IkBn*_NFkBn) * dt;
	dA20 = (_PARAM(PARAM_GR_c4)*_A20t-_PARAM(PARAM_GR_c5)*_A20) * dt;
	dA20t = (_PARAM(PARAM_GR_c1)*_GA20-_PARAM(PARAM_GR_c3)*_A20t) * dt;
	dIkB = (-_PARAM(PARAM_GR_a2)*_IKKa*_IkB-_PARAM(PARAM_GR_a1)*_IkB*_NFkBc+_PARAM(PARAM_GR_c4)*_IkBt-_PARAM(PARAM_GR_c5a)*_IkB-_PARAM(PARAM_GR_i1a)*_IkB+_PARAM(PARAM_GR_e1a)*_IkBn) * dt;
	dIkBn = (-_PARAM(PARAM_GR_a1)*kv*_IkBn*_NFkBn+_PARAM(PARAM_GR_i1a)*_IkB-_PARAM(PARAM_GR_e1a)*_IkBn) * dt;
	dIkBt = (_PARAM(PARAM_GR_c1)*_GIkB-_PARAM(PARAM_GR_c3)*_IkBt) * dt;
	dNFkB_IkB = (_PARAM(PARAM_GR_a1)*_IkB*_NFkBc-_PARAM(PARAM_GR_c6a)*_NFkB_IkB-_PARAM(PARAM_GR_a3)*_IKKa*_NFkB_IkB+_PARAM(PARAM_GR_e2a)*_NFkB_IkBn) * dt;
	dNFkB_IkBn = (_PARAM(PARAM_GR_a1)*kv*_IkBn*_NFkBn-_PARAM(PARAM_GR_e2a)*_NFkB_IkBn) * dt;
	dGA20 = (_PARAM(PARAM_GR_q1)*_NFkBn*(2-_GA20)-_PARAM(PARAM_GR_q2)*_IkBn*_GA20) * dt;
	dGIkB = (_PARAM(PARAM_GR_q1)*_NFkBn*(2-_GIkB)-_PARAM(PARAM_GR_q2)*_IkBn*_GIkB) * dt;
	dGR = (_PARAM(PARAM_GR_q1r)*_NFkBn*(2-_GR)-(_PARAM(PARAM_GR_q2rr)+_PARAM(PARAM_GR_q2r)*_IkBn)*_GR) * dt;
	dchemt = (_c1rrChemTNF+_c1rChem*_GR-_PARAM(PARAM_GR_c3rChem)*_chemt) * dt;
	dchem = (_PARAM(PARAM_GR_c4Chem)*_chemt-_PARAM(PARAM_GR_c5Chem)*_chem-_PARAM(PARAM_GR_e3Chem)*_chem) * dt;
	dTNFt = (_c1rrChemTNF+_c1rTNF*_GR-_PARAM(PARAM_GR_c3rTNF)*_TNFt) * dt;
	dTNF = (_PARAM(PARAM_GR_c4TNF)*_TNFt-_PARAM(PARAM_GR_c5TNF)*_TNF-_PARAM(PARAM_GR_e3TNF)*_TNF) * dt;
	dACTt = (_PARAM(PARAM_GR_c1rrACT)+_PARAM(PARAM_GR_c1r)*_GR-_PARAM(PARAM_GR_c3rACT)*_ACTt) * dt;
	dACT = (_PARAM(PARAM_GR_c4ACT)*_ACTt-_PARAM(PARAM_GR_c5ACT)*_ACT) * dt;
	dIAPt = (_PARAM(PARAM_GR_c1rrIAP)+_PARAM(PARAM_GR_c1r)*_GR-_PARAM(PARAM_GR_c3rIAP)*_IAPt) * dt;
	dIAP = (_PARAM(PARAM_GR_c4IAP)*_IAPt-_PARAM(PARAM_GR_c5IAP)*_IAP) * dt;
	
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
	
	// secrete chemokines
	cell.incCCL2(_PARAM(PARAM_GR_e3Chem) * _chem * dt);
	cell.incCCL5(_PARAM(PARAM_GR_e3Chem) * _chem * dt);
	cell.incCXCL9(2 * _PARAM(PARAM_GR_e3Chem) * _chem * dt);
	
	_IKKKa += dIKKKa;
	_IKKn += dIKKn; 
	_IKKa += dIKKa; 
	_IKKi += dIKKi; 
	_IkBp += dIkBp; 
	_NFkB_IkBp += dNFkB_IkBp; 
	_NFkBc += dNFkBc; 
	_NFkBn += dNFkBn; 
	_A20 += dA20;
	_A20t += dA20t; 
	_IkB += dIkB;
	_IkBn += dIkBn; 
	_IkBt += dIkBt; 
	_NFkB_IkB += dNFkB_IkB;
	_NFkB_IkBn += dNFkB_IkBn;
	_GA20 += dGA20; 
	_GIkB += dGIkB;  
	_GR += dGR; 
	_chemt += dchemt; 
	_chem += dchem;
	_TNFt += dTNFt; 
	_TNF += dTNF;
	_ACTt += dACTt; 
	_ACT += dACT;
	_normalizedACT = _ACT*_PARAM(PARAM_GR_c3rACT);
	_IAPt += dIAPt; 
	_IAP += dIAP;
	_normalizedIAP = _IAP*_PARAM(PARAM_GR_c3rIAP);
	
	if (_mTNF < 0 || _surfTNFR1 < 0 || _surfBoundTNFR1 < 0 || _surfTNFR2 < 0 || _surfBoundTNFR2 < 0)
		std::cout << "Error: Negative Value of Species in TNF/TNFR dynamics" << std::endl;
	if (_IKKKa < 0 || _IKKn < 0 || _IKKa < 0 || _IKKa < 0 || _IKKi < 0 || _IkBp < 0 || _NFkB_IkBp < 0 || 
		_NFkBc < 0 || _NFkBn < 0 || _A20 < 0 || _A20t < 0 || _IkB < 0 || _IkBn < 0 || _IkBt < 0 ||
		_NFkB_IkB < 0 || _NFkB_IkBn < 0 || _GA20 < 0 || _GIkB < 0 || _GR < 0 || _chemt < 0 || _chem < 0 ||
		_TNFt < 0 || _TNF < 0 || _ACTt < 0 || _ACT < 0 || _IAPt < 0 || _IAP < 0)
		std::cout << "Error: Negative Value of Species in NFkB dynamics" << std::endl;
    
    //cout << "Debug: Running TNF and NFkB dynamics" << std::endl;
}



void Mac::solveTNFandIL10(GrGrid& grid, GrStat& stats, double dt, double currenttime)
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
        IkmRNA = _kmRNA * ((_kSynth/_kmRNA) + ((1.0 - (_kSynth/_kmRNA))/(1.0 + pow(2.7183, ((_surfBoundIL10R - _PARAM(PARAM_GR_LINK_RNA_GAMMA))/_PARAM(PARAM_GR_LINK_RNA_DELTA))))));
		
		if (IkmRNA <= 0.75)
		{
			stats.incNrOfCellsInhibited();
		}
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
    
   // cout << dsurfBoundIL10R << std::endl;
    
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
	
	if (currenttime == 99)
	{
		if (_state == MAC_INFECTED && _NFkB > 0 && _deactivationTime == -1)
		{
			cout << "NFkB Infected    " << IkmRNA << "     " << "Bound IL10    "<< _surfBoundIL10R << std::endl;
		}
		
		if (_state == MAC_INFECTED && _NFkB < 1 && _deactivationTime == -1)
		{
			cout << "Infected    " << IkmRNA << "     " << "Bound IL10    "<< _surfBoundIL10R << std::endl;
		}
		
		if (_state == MAC_ACTIVE && _NFkB > 0 && _deactivationTime == -1)
		{
			cout << "NFkB Active    " << IkmRNA << "     " << "Bound IL10    "<< _surfBoundIL10R << std::endl;
		}
		
		if (_state == MAC_ACTIVE && _NFkB < 1 && _deactivationTime == -1)
		{
			cout << "Active    " << IkmRNA << "     " << "Bound IL10    "<< _surfBoundIL10R << std::endl;
		}
	
	}
	cell.setTNF(Nav * vol * tnf);
	cell.setShedTNFR2(Nav * vol * shedtnfr2);
    cell.setIL10(Nav * vol * il10);
	
    
	if (_mTNF < 0 || _surfTNFR1 < 0 || _surfBoundTNFR1 < 0 || _surfTNFR2 < 0 || _surfBoundTNFR2 < 0 || _mTNFRNA < 0)
		std::cout << "Error: Negative Value of Species in TNF/TNFR dynamics" << std::endl;
    
    if (_surfIL10R < 0 || _surfBoundIL10R < 0)
        std::cout << "Error: Negative value of species in IL10/IL10R dynamics" << std::endl;
    
    //cout << "Debug: Running TNF and IL10 dynamics" << std::endl;
    
}


void Mac::solveTNFandIL10andNFkB(GrGrid& grid, double dt)
{
   	GridCell& cell = grid(_row, _col);
	
	double koff1 = _PARAM(PARAM_GR_K_ON1) * _PARAM(PARAM_GR_KD1);
	double koff2 = _PARAM(PARAM_GR_K_ON2) * _PARAM(PARAM_GR_KD2);
	double kdeg = 4.58e-4;
	double density = 1.25e11; // used for conversion of conc. unit (M -> #/cell) based on cell and microcompartment volumes 
	double Nav = 6.02e23; // Avogadro Number
	double vol = 8.0e-12; // volume of a cell in liter
	double kv = 5; //ratio of cytoplasmic to nuclear volume
	
	double tnf = cell.getTNF() / (Nav * vol);
	double shedtnfr2 = cell.getShedTNFR2() / (Nav * vol);
    double il10 = cell.getIL10() / (Nav * vol);
	
	double dmTNF;
	double dsurfTNFR1;
	double dsurfTNFR2;
	double dsurfBoundTNFR1;
	double dsurfBoundTNFR2;
	double dintBoundTNFR1;
	double dintBoundTNFR2;
	double dsTNF;
	double dshedTNFR2;
	
	double dIKKKa; 
	double dIKKn; 
	double dIKKa; 
	double dIKKi; 
	double dIkBp; 
	double dNFkB_IkBp; 
	double dNFkBc; 
	double dNFkBn; 
	double dA20;
	double dA20t; 
	double dIkB;
	double dIkBn; 
	double dIkBt; 
	double dNFkB_IkB;
	double dNFkB_IkBn;
	double dGA20;
	double dGIkB;
	double dGR; 
	double dchemt; 
	double dchem; 
	double dTNFt;
	double dTNF; 
	double dACTt; 
	double dACT; 
	double dIAPt; 
	double dIAP;
    
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

    
     // CURRENTLY NOT LINKED WITH IL10 SINCE NFKB DYNAMICS HAVE NOT BEEN LOOKED AT
    
	dmTNF = (_PARAM(PARAM_GR_e3TNF)*_TNF - _kTACE * _mTNF) * dt;
	dsurfTNFR1 = (_vTNFR1 - _PARAM(PARAM_GR_K_ON1) * tnf * _surfTNFR1 + (koff1+kdeg) * _surfBoundTNFR1 - _PARAM(PARAM_GR_K_T1) * _surfTNFR1 + _PARAM(PARAM_GR_K_REC1) * _intBoundTNFR1) * dt;
	dsurfTNFR2 = (_vTNFR2 - _PARAM(PARAM_GR_K_ON2) * tnf * _surfTNFR2 + (koff2+kdeg) * _surfBoundTNFR2 - _PARAM(PARAM_GR_K_T2) * _surfTNFR2 + _PARAM(PARAM_GR_K_REC2) * _intBoundTNFR2) * dt;
	dsurfBoundTNFR1 = (_PARAM(PARAM_GR_K_ON1) * tnf * _surfTNFR1 - (koff1+kdeg) * _surfBoundTNFR1 - _PARAM(PARAM_GR_K_INT1) * _surfBoundTNFR1) * dt;
	dsurfBoundTNFR2 = (_PARAM(PARAM_GR_K_ON2) * tnf * _surfTNFR2 - (koff2+kdeg) * _surfBoundTNFR2 - _PARAM(PARAM_GR_K_INT2) * _surfBoundTNFR2 - _PARAM(PARAM_GR_K_SHED) * _surfBoundTNFR2) * dt;
	dintBoundTNFR1 = (_PARAM(PARAM_GR_K_INT1) * _surfBoundTNFR1 - _PARAM(PARAM_GR_K_DEG1) * _intBoundTNFR1 - _PARAM(PARAM_GR_K_REC1) * _intBoundTNFR1) * dt;
	dintBoundTNFR2 = (_PARAM(PARAM_GR_K_INT2) * _surfBoundTNFR2 - _PARAM(PARAM_GR_K_DEG2) * _intBoundTNFR2 - _PARAM(PARAM_GR_K_REC2) * _intBoundTNFR2) * dt;
	dsTNF = ((density/Nav) * (_kTACE * _mTNF - _PARAM(PARAM_GR_K_ON1) * tnf * _surfTNFR1 + koff1 * _surfBoundTNFR1 - _PARAM(PARAM_GR_K_ON2) * tnf * _surfTNFR2 + koff2 * _surfBoundTNFR2)) * dt; 
	dshedTNFR2 = ((density/Nav) * _PARAM(PARAM_GR_K_SHED) * _surfBoundTNFR2) * dt;
	
	// NF-kB dynamics model equations
	dIKKKa = (_PARAM(PARAM_GR_ka)*_surfBoundTNFR1*(_PARAM(PARAM_GR_KN)-_IKKKa)*_PARAM(PARAM_GR_kA20)/(_PARAM(PARAM_GR_kA20)+_A20)-_PARAM(PARAM_GR_ki)*_IKKKa) * dt;
	dIKKn = (_PARAM(PARAM_GR_k4)*(_PARAM(PARAM_GR_KNN)-_IKKn-_IKKa-_IKKi)-_PARAM(PARAM_GR_k1)*pow(_IKKKa,2.0)*_IKKn) * dt;
	dIKKa = (_PARAM(PARAM_GR_k1)*pow(_IKKKa,2.0)*_IKKn-_PARAM(PARAM_GR_k3)*_IKKa*(_PARAM(PARAM_GR_k2)+_A20)/_PARAM(PARAM_GR_k2)) * dt;
	dIKKi = (_PARAM(PARAM_GR_k3)*_IKKa*(_PARAM(PARAM_GR_k2)+_A20)/_PARAM(PARAM_GR_k2)-_PARAM(PARAM_GR_k4)*_IKKi) * dt;
	dIkBp = (_PARAM(PARAM_GR_a2)*_IKKa*_IkB-_PARAM(PARAM_GR_tp)*_IkBp) * dt;
	dNFkB_IkBp = (_PARAM(PARAM_GR_a3)*_IKKa*_NFkB_IkB-_PARAM(PARAM_GR_tp)*_NFkB_IkBp) * dt;
	dNFkBc = (_PARAM(PARAM_GR_c6a)*_NFkB_IkB-_PARAM(PARAM_GR_a1)*_NFkBc*_IkB+_PARAM(PARAM_GR_tp)*_NFkB_IkBp-_PARAM(PARAM_GR_i1)*_NFkBc) * dt;
	dNFkBn = (_PARAM(PARAM_GR_i1)*_NFkBc-_PARAM(PARAM_GR_a1)*kv*_IkBn*_NFkBn) * dt;
	dA20 = (_PARAM(PARAM_GR_c4)*_A20t-_PARAM(PARAM_GR_c5)*_A20) * dt;
	dA20t = (_PARAM(PARAM_GR_c1)*_GA20-_PARAM(PARAM_GR_c3)*_A20t) * dt;
	dIkB = (-_PARAM(PARAM_GR_a2)*_IKKa*_IkB-_PARAM(PARAM_GR_a1)*_IkB*_NFkBc+_PARAM(PARAM_GR_c4)*_IkBt-_PARAM(PARAM_GR_c5a)*_IkB-_PARAM(PARAM_GR_i1a)*_IkB+_PARAM(PARAM_GR_e1a)*_IkBn) * dt;
	dIkBn = (-_PARAM(PARAM_GR_a1)*kv*_IkBn*_NFkBn+_PARAM(PARAM_GR_i1a)*_IkB-_PARAM(PARAM_GR_e1a)*_IkBn) * dt;
	dIkBt = (_PARAM(PARAM_GR_c1)*_GIkB-_PARAM(PARAM_GR_c3)*_IkBt) * dt;
	dNFkB_IkB = (_PARAM(PARAM_GR_a1)*_IkB*_NFkBc-_PARAM(PARAM_GR_c6a)*_NFkB_IkB-_PARAM(PARAM_GR_a3)*_IKKa*_NFkB_IkB+_PARAM(PARAM_GR_e2a)*_NFkB_IkBn) * dt;
	dNFkB_IkBn = (_PARAM(PARAM_GR_a1)*kv*_IkBn*_NFkBn-_PARAM(PARAM_GR_e2a)*_NFkB_IkBn) * dt;
	dGA20 = (_PARAM(PARAM_GR_q1)*_NFkBn*(2-_GA20)-_PARAM(PARAM_GR_q2)*_IkBn*_GA20) * dt;
	dGIkB = (_PARAM(PARAM_GR_q1)*_NFkBn*(2-_GIkB)-_PARAM(PARAM_GR_q2)*_IkBn*_GIkB) * dt;
	dGR = (_PARAM(PARAM_GR_q1r)*_NFkBn*(2-_GR)-(_PARAM(PARAM_GR_q2rr)+_PARAM(PARAM_GR_q2r)*_IkBn)*_GR) * dt;
	dchemt = (_c1rrChemTNF+_c1rChem*_GR-_PARAM(PARAM_GR_c3rChem)*_chemt) * dt;
	dchem = (_PARAM(PARAM_GR_c4Chem)*_chemt-_PARAM(PARAM_GR_c5Chem)*_chem-_PARAM(PARAM_GR_e3Chem)*_chem) * dt;
	dTNFt = (_c1rrChemTNF+_c1rTNF*_GR-_PARAM(PARAM_GR_c3rTNF)*_TNFt) * dt;
	dTNF = (_PARAM(PARAM_GR_c4TNF)*_TNFt-_PARAM(PARAM_GR_c5TNF)*_TNF-_PARAM(PARAM_GR_e3TNF)*_TNF) * dt;
	dACTt = (_PARAM(PARAM_GR_c1rrACT)+_PARAM(PARAM_GR_c1r)*_GR-_PARAM(PARAM_GR_c3rACT)*_ACTt) * dt;
	dACT = (_PARAM(PARAM_GR_c4ACT)*_ACTt-_PARAM(PARAM_GR_c5ACT)*_ACT) * dt;
	dIAPt = (_PARAM(PARAM_GR_c1rrIAP)+_PARAM(PARAM_GR_c1r)*_GR-_PARAM(PARAM_GR_c3rIAP)*_IAPt) * dt;
	dIAP = (_PARAM(PARAM_GR_c4IAP)*_IAPt-_PARAM(PARAM_GR_c5IAP)*_IAP) * dt;
	
    // IL10 differential equations
    dsIL10 = (density/Nav) * _kISynth + ((density/Nav) * (_PARAM(PARAM_GR_I_K_OFF) * _surfBoundIL10R - _PARAM(PARAM_GR_I_K_ON) * _surfIL10R * il10)) * dt;
	dsurfIL10R = (_vIL10R - _PARAM(PARAM_GR_I_K_ON) * _surfIL10R * il10 + _PARAM(PARAM_GR_I_K_OFF) * _surfBoundIL10R - _PARAM(PARAM_GR_I_K_T) * _surfIL10R) * dt;
	dsurfBoundIL10R = (_PARAM(PARAM_GR_I_K_ON) * _surfIL10R * il10 - _PARAM(PARAM_GR_I_K_OFF) * _surfBoundIL10R - _PARAM(PARAM_GR_I_K_INT) * _surfBoundIL10R) * dt;
    // end of IL10 differential equations
    
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
    
	// secrete chemokines
	cell.incCCL2(_PARAM(PARAM_GR_e3Chem) * _chem * dt);
	cell.incCCL5(_PARAM(PARAM_GR_e3Chem) * _chem * dt);
	cell.incCXCL9(2 * _PARAM(PARAM_GR_e3Chem) * _chem * dt);
	
	_IKKKa += dIKKKa;
	_IKKn += dIKKn; 
	_IKKa += dIKKa; 
	_IKKi += dIKKi; 
	_IkBp += dIkBp; 
	_NFkB_IkBp += dNFkB_IkBp; 
	_NFkBc += dNFkBc; 
	_NFkBn += dNFkBn; 
	_A20 += dA20;
	_A20t += dA20t; 
	_IkB += dIkB;
	_IkBn += dIkBn; 
	_IkBt += dIkBt; 
	_NFkB_IkB += dNFkB_IkB;
	_NFkB_IkBn += dNFkB_IkBn;
	_GA20 += dGA20; 
	_GIkB += dGIkB;  
	_GR += dGR; 
	_chemt += dchemt; 
	_chem += dchem;
	_TNFt += dTNFt; 
	_TNF += dTNF;
	_ACTt += dACTt; 
	_ACT += dACT;
	_normalizedACT = _ACT*_PARAM(PARAM_GR_c3rACT);
	_IAPt += dIAPt; 
	_IAP += dIAP;
	_normalizedIAP = _IAP*_PARAM(PARAM_GR_c3rIAP);
	
	if (_mTNF < 0 || _surfTNFR1 < 0 || _surfBoundTNFR1 < 0 || _surfTNFR2 < 0 || _surfBoundTNFR2 < 0)
		std::cout << "Error: Negative Value of Species in TNF/TNFR dynamics" << std::endl;
    if (_surfIL10R < 0 || _surfBoundIL10R < 0)
        std::cout << "Error: Negative value of species in IL10/IL10R dynamics" << std::endl;
    if (_IKKKa < 0 || _IKKn < 0 || _IKKa < 0 || _IKKa < 0 || _IKKi < 0 || _IkBp < 0 || _NFkB_IkBp < 0 || 
		_NFkBc < 0 || _NFkBn < 0 || _A20 < 0 || _A20t < 0 || _IkB < 0 || _IkBn < 0 || _IkBt < 0 ||
		_NFkB_IkB < 0 || _NFkB_IkBn < 0 || _GA20 < 0 || _GIkB < 0 || _GR < 0 || _chemt < 0 || _chem < 0 ||
		_TNFt < 0 || _TNF < 0 || _ACTt < 0 || _ACT < 0 || _IAPt < 0 || _IAP < 0)
		std::cout << "Error: Negative Value of Species in NFkB dynamics" << std::endl;
    
     //cout << "Debug: Running TNF and IL10 and NFkB dynamics" << std::endl;
}


void Mac::solveIL10(GrGrid& grid, double dt)
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
    
     //cout << "Debug: Running IL10 dynamics" << std::endl;
    
}

void Mac::solveDegradation(GrGrid& grid, double dt, bool tnfrDynamics, bool il10rDynamics)
{
    GridCell& cell = grid(_row, _col);
    
    double Nav = 6.02e23; // Avagadros #
    double vol = 8.0e-12; // volume of a cell in L
    
    if (!tnfrDynamics) {
        
        // simulate the effect of TNF internalization by cells in the form of degradation. Only for TNF
        double dtnf;
        double tnf = cell.getTNF();

        dtnf = -_PARAM(PARAM_GR_K_INT1) * (tnf / (tnf + _PARAM(PARAM_GR_KD1) * Nav * vol)) * _PARAM(PARAM_GR_MEAN_TNFR1_MAC) * dt * 0.4;
        tnf += dtnf;
        
        cell.setTNF(tnf);
    }
    
    if (!il10rDynamics) {
        
        double dil10;
        double il10 = cell.getIL10();
        
        // simulate the effect of IL10 internalization in the form of degradation. Only for IL10
        
        dil10 = -_PARAM(PARAM_GR_I_K_INT) * (il10 / (il10 + _PARAM(PARAM_GR_I_KD) * Nav * vol)) * _PARAM(PARAM_GR_I_IL10R_MAC) * dt * _PARAM(PARAM_GR_I_MOD);
        il10 += dil10;
        
        cell.setIL10(il10);
    }
    
}


void Mac::solveNFkBODEsEquilibrium(double dt)
{
	double kv = 5; //ratio of cytoplasmic to nuclear volume
	
	double dIKKKa; 
	double dIKKn; 
	double dIKKa; 
	double dIKKi; 
	double dIkBp; 
	double dNFkB_IkBp; 
	double dNFkBc; 
	double dNFkBn; 
	double dA20;
	double dA20t; 
	double dIkB;
	double dIkBn; 
	double dIkBt; 
	double dNFkB_IkB;
	double dNFkB_IkBn;
	double dGA20;
	double dGIkB;
	double dGR; 
	double dchemt; 
	double dchem; 
	double dTNFt;
	double dTNF; 
	double dACTt; 
	double dACT; 
	double dIAPt; 
	double dIAP;
	
	// NF-kB dynamics model equations
	//dIKKKa = (_PARAM(PARAM_GR_ka)*_surfBoundTNFR1*(_PARAM(PARAM_GR_KN)-_IKKKa)*_PARAM(PARAM_GR_kA20)/(_PARAM(PARAM_GR_kA20)+_A20)-_PARAM(PARAM_GR_ki)*_IKKKa) * dt;
	dIKKKa = 0.0;
	//dIKKn = (_PARAM(PARAM_GR_k4)*(_PARAM(PARAM_GR_KNN)-_IKKn-_IKKa-_IKKi)-_PARAM(PARAM_GR_k1)*pow(_IKKKa,2.0)*_IKKn) * dt;
	dIKKn = 0.0;
	//dIKKa = (_PARAM(PARAM_GR_k1)*pow(_IKKKa,2.0)*_IKKn-_PARAM(PARAM_GR_k3)*_IKKa*(_PARAM(PARAM_GR_k2)+_A20)/_PARAM(PARAM_GR_k2)) * dt;
	dIKKa = 0.0;
	//dIKKi = (_PARAM(PARAM_GR_k3)*_IKKa*(_PARAM(PARAM_GR_k2)+_A20)/_PARAM(PARAM_GR_k2)-_PARAM(PARAM_GR_k4)*_IKKi) * dt;
	dIKKi = 0.0;
	//dIkBp = (_PARAM(PARAM_GR_a2)*_IKKa*_IkB-_PARAM(PARAM_GR_tp)*_IkBp) * dt;
	dIkBp = 0.0;
	//dNFkB_IkBp = (_PARAM(PARAM_GR_a3)*_IKKa*_NFkB_IkB-_PARAM(PARAM_GR_tp)*_NFkB_IkBp) * dt;
	dNFkB_IkBp = 0.0;
	dNFkBc = (_PARAM(PARAM_GR_c6a)*_NFkB_IkB-_PARAM(PARAM_GR_a1)*_NFkBc*_IkB+_PARAM(PARAM_GR_tp)*_NFkB_IkBp-_PARAM(PARAM_GR_i1)*_NFkBc) * dt;
	dNFkBn = (_PARAM(PARAM_GR_i1)*_NFkBc-_PARAM(PARAM_GR_a1)*kv*_IkBn*_NFkBn) * dt;
	dA20 = (_PARAM(PARAM_GR_c4)*_A20t-_PARAM(PARAM_GR_c5)*_A20) * dt;
	dA20t = (_PARAM(PARAM_GR_c1)*_GA20-_PARAM(PARAM_GR_c3)*_A20t) * dt;
	dIkB = (-_PARAM(PARAM_GR_a2)*_IKKa*_IkB-_PARAM(PARAM_GR_a1)*_IkB*_NFkBc+_PARAM(PARAM_GR_c4)*_IkBt-_PARAM(PARAM_GR_c5a)*_IkB-_PARAM(PARAM_GR_i1a)*_IkB+_PARAM(PARAM_GR_e1a)*_IkBn) * dt;
	dIkBn = (-_PARAM(PARAM_GR_a1)*kv*_IkBn*_NFkBn+_PARAM(PARAM_GR_i1a)*_IkB-_PARAM(PARAM_GR_e1a)*_IkBn) * dt;
	dIkBt = (_PARAM(PARAM_GR_c1)*_GIkB-_PARAM(PARAM_GR_c3)*_IkBt) * dt;
	dNFkB_IkB = (_PARAM(PARAM_GR_a1)*_IkB*_NFkBc-_PARAM(PARAM_GR_c6a)*_NFkB_IkB-_PARAM(PARAM_GR_a3)*_IKKa*_NFkB_IkB+_PARAM(PARAM_GR_e2a)*_NFkB_IkBn) * dt;
	dNFkB_IkBn = (_PARAM(PARAM_GR_a1)*kv*_IkBn*_NFkBn-_PARAM(PARAM_GR_e2a)*_NFkB_IkBn) * dt;
	dGA20 = (_PARAM(PARAM_GR_q1)*_NFkBn*(2-_GA20)-_PARAM(PARAM_GR_q2)*_IkBn*_GA20) * dt;
	dGIkB = (_PARAM(PARAM_GR_q1)*_NFkBn*(2-_GIkB)-_PARAM(PARAM_GR_q2)*_IkBn*_GIkB) * dt;
	dGR = (_PARAM(PARAM_GR_q1r)*_NFkBn*(2-_GR)-(_PARAM(PARAM_GR_q2rr)+_PARAM(PARAM_GR_q2r)*_IkBn)*_GR) * dt;
	//dchemt = (_c1rrChemTNF+_c1rChem*_GR-_PARAM(PARAM_GR_c3rChem)*_chemt) * dt;
	dchemt = 0.0;
	//dchem = (_PARAM(PARAM_GR_c4Chem)*_chemt-_PARAM(PARAM_GR_c5Chem)*_chem-_PARAM(PARAM_GR_e3Chem)*_chem) * dt; 
	dchem = 0.0;
	//dTNFt = (_c1rrChemTNF+_c1rTNF*_GR-_PARAM(PARAM_GR_c3rTNF)*_TNFt) * dt;
	dTNFt = 0.0;
	//dTNF = (_PARAM(PARAM_GR_c4TNF)*_TNFt-_PARAM(PARAM_GR_c5TNF)*_TNF-_PARAM(PARAM_GR_e3TNF)*_TNF) * dt;
	dTNF = 0.0;
	//dACTt = (_PARAM(PARAM_GR_c1rrACT)+_PARAM(PARAM_GR_c1r)*_GR-_PARAM(PARAM_GR_c3rACT)*_ACTt) * dt;
	dACTt = 0.0;
	//dACT = (_PARAM(PARAM_GR_c4ACT)*_ACTt-_PARAM(PARAM_GR_c5ACT)*_ACT) * dt;
	dACT = 0.0;
	//dIAPt = (_PARAM(PARAM_GR_c1rrIAP)+_PARAM(PARAM_GR_c1r)*_GR-_PARAM(PARAM_GR_c3rIAP)*_IAPt) * dt;
	dIAPt = 0.0;
	//dIAP = (_PARAM(PARAM_GR_c4IAP)*_IAPt-_PARAM(PARAM_GR_c5IAP)*_IAP) * dt;
	dIAP = 0.0;
	
	/*
	 some species as shown above, during the equilibrium stage in the absence of TNF, won't change from their initial 
	 values, so their differential equations are set to zero.
	 */
	
	_IKKKa += dIKKKa;
	_IKKn += dIKKn; 
	_IKKa += dIKKa; 
	_IKKi += dIKKi; 
	_IkBp += dIkBp; 
	_NFkB_IkBp += dNFkB_IkBp; 
	_NFkBc += dNFkBc; 
	_NFkBn += dNFkBn; 
	_A20 += dA20;
	_A20t += dA20t; 
	_IkB += dIkB;
	_IkBn += dIkBn; 
	_IkBt += dIkBt; 
	_NFkB_IkB += dNFkB_IkB;
	_NFkB_IkBn += dNFkB_IkBn;
	_GA20 += dGA20; 
	_GIkB += dGIkB;  
	_GR += dGR; 
	_chemt += dchemt; 
	_chem += dchem;
	_TNFt += dTNFt; 
	_TNF += dTNF;
	_ACTt += dACTt; 
	_ACT += dACT;
	_normalizedACT = _ACT*_PARAM(PARAM_GR_c3rACT);
	_IAPt += dIAPt; 
	_IAP += dIAP;
	_normalizedIAP = _IAP*_PARAM(PARAM_GR_c3rIAP);
	/*	
	 if (_IKKKa < 0 || _IKKn < 0 || _IKKa < 0 || _IKKa < 0 || _IKKi < 0 || _IkBp < 0 || _NFkB_IkBp < 0 || 
	 _NFkBc < 0 || _NFkBn < 0 || _A20 < 0 || _A20t < 0 || _IkB < 0 || _IkBn < 0 || _IkBt < 0 ||
	 _NFkB_IkB < 0 || _NFkB_IkBn < 0 || _GA20 < 0 || _GIkB < 0 || _GR < 0 || _chemt < 0 || _chem < 0 ||
	 _TNFt < 0 || _TNF < 0 || _ACTt < 0 || _ACT < 0 || _IAPt < 0 || _IAP < 0)
	 std::cout << "Error: Negative Value of Species in NFkB dynamics" << std::endl;
	 */
}

void Mac::handleResting(const int time, GrGrid& grid, GrStat& stats, bool nfkbDynamics)
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
	if (nfkbDynamics) // macrophage activation dynamics follow NFkB dynamics
	{
		if (_stat1)
		{
			if (getExtMtbInMoore(grid) > _PARAM(PARAM_MAC_THRESHOLD_NFKB_EXTMTB))
			{
				_intMtb = 0;
				_activationTime = time;
				_deathTime = time + _PARAM(PARAM_MAC_A_AGE);
				_nextState = MAC_ACTIVE;
			}
			else if (_normalizedACT > _PARAM(PARAM_GR_ACT_THRESHOLD) &&
					 g_Rand.getReal() < 1 - pow(2.7183, -_PARAM(PARAM_GR_ACT_K) * (_normalizedACT - _PARAM(PARAM_GR_ACT_THRESHOLD))))
			{
				_intMtb = 0;
				_activationTime = time;
				_deathTime = time + _PARAM(PARAM_MAC_A_AGE);
				_nextState = MAC_ACTIVE;
				stats.incRestingMacActivationTNF();
			}
		}
	}
	else if (_stat1 && _NFkB)
	{
		_intMtb = 0;
		_activationTime = time;
		_deathTime = time + _PARAM(PARAM_MAC_A_AGE);
		_nextState = MAC_ACTIVE;
		stats.incRestingMacActivationTNF();
	}
}

void Mac::handleInfected(const int time, GrGrid& grid, GrStat& stats, bool nfkbDynamics)
{
	GridCell& cell = grid(_row, _col);

	// intracellular bacteria reproduce
	_intMtb *= getIntMtbGrowthRate(time);

	// uptake extracellular bacteria
	// The probability to uptake decreases linearly as the number of intra-cellular bacteria increases to the threshold
	// for becoming chronically infected. At that threshold the probability becomes 0.
	// The probability of uptake is the compliment of the probability of killing.
	double baseProbExtMtbUptake = (1.0 - _PARAM(PARAM_MAC_PROB_KILL_R_EXTMTB))/2.0;
	double probExtMtbUptake = (baseProbExtMtbUptake *  (1.0 - (_intMtb / _PARAM(PARAM_MAC_THRESHOLD_BECOME_CI_INTMTB)) ));
	if (cell.getExtMtb() > 0 && g_Rand.getReal() < probExtMtbUptake)
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
					_nextState = MAC_ACTIVE;
				}
				else if (_normalizedACT > _PARAM(PARAM_GR_ACT_THRESHOLD) &&
						 g_Rand.getReal() < 1 - pow(2.7183, -_PARAM(PARAM_GR_ACT_K) * (_normalizedACT - _PARAM(PARAM_GR_ACT_THRESHOLD))))
				{
					_intMtb = 0;
					_activationTime = time;
					_deathTime = time + _PARAM(PARAM_MAC_A_AGE);
					_nextState = MAC_ACTIVE;
					stats.incInfMacActivationTNF();
				}
				else 
				{
					_nextState = MAC_INFECTED;
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
				_nextState = MAC_ACTIVE;
				stats.incInfMacActivationTNF();
			}
			else
			{
				_nextState = MAC_INFECTED;
			}
		}
	}
}

void Mac::handleChronicallyInfected(const int time, GrGrid& grid, GrStat&)
{
	GridCell& cell = grid(_row, _col);

	// intracellular bacteria reproduce
	_intMtb *= getIntMtbGrowthRate(time);

	if (_intMtb >= _PARAM(PARAM_MAC_THRESHOLD_BURST_CI_INTMTB))
	{
		// burst, all intracellular bacteria disperse to the Moore neighborhood
		disperseMtb(grid, 1.0);
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
	
	out << _IKKKa << std::endl;
	out << _IKKn << std::endl;
	out << _IKKa << std::endl;
	out << _IKKi << std::endl;
	out << _IkBp << std::endl;
	out << _NFkB_IkBp << std::endl;
	out << _NFkBc << std::endl;
	out << _NFkBn << std::endl;
	out << _A20 << std::endl;
	out << _A20t << std::endl;
	out << _IkB << std::endl;
	out << _IkBn << std::endl;
	out << _IkBt << std::endl;
	out << _NFkB_IkB << std::endl;
	out << _NFkB_IkBn << std::endl;
	out << _GA20 << std::endl;
	out << _GIkB << std::endl;
	out << _GR << std::endl;
	out << _c1rrChemTNF << std::endl;
	out << _c1rChem << std::endl;
	out << _c1rTNF << std::endl;
	out << _chemt << std::endl;
	out << _chem << std::endl;
	out << _TNFt << std::endl;
	out << _TNF << std::endl;
	out << _ACTt << std::endl;
	out << _ACT << std::endl;
	out << _normalizedACT << std::endl;
	out << _IAPt << std::endl;
	out << _IAP << std::endl;
	out << _normalizedIAP << std::endl;

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

	in >> intVal;
	_state = (MacState) intVal;

	in >> intVal;
	_nextState = (MacState) intVal;

	in >> _intMtb;
	in >> _NFkB;
	in >> _stat1;
	in >> _ICOS;
	in >> _activationTime;
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
    
	in >> _IKKKa;
	in >> _IKKn;
	in >> _IKKa;
	in >> _IKKi;
	in >> _IkBp;
	in >> _NFkB_IkBp;
	in >> _NFkBc;
	in >> _NFkBn;
	in >> _A20;
	in >> _A20t;
	in >> _IkB;
	in >> _IkBn;
	in >> _IkBt;
	in >> _NFkB_IkB;
	in >> _NFkB_IkBn;
	in >> _GA20;
	in >> _GIkB;
	in >> _GR;
	in >> _c1rrChemTNF;
	in >> _c1rChem;
	in >> _c1rTNF;
	in >> _chemt;
	in >> _chem;
	in >> _TNFt;
	in >> _TNF;
	in >> _ACTt;
	in >> _ACT;
	in >> _normalizedACT;
	in >> _IAPt;
	in >> _IAP;
	in >> _normalizedIAP;
	
	if (!Serialization::readFooter(in, Mac::_ClassName))
	{
		exit(1);
	}
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

	// count the number of extMtb cells in the Moore neighborhood
	// that are in the specified state
	for (int i = -1; i <= 1; i++)
		for (int j = -1; j <= 1; j++)
			mtb += grid(MOD_ROW(_row + i), MOD_COL(_col + j)).getExtMtb();

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
    	Pos pos = compartmentOrdinalToCoordinates(permutation[i]);
    	int row = pos.first;
    	int col = pos.second;

    	GridCell& cell = grid(row, col);
    	cell.incExtMtb(dMtb);
    }

    _intMtb -= mtb;
}

