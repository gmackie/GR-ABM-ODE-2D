/*
 * agent.cpp
 *
 *  Created on: 02-nov-2009
 *      Author: M. El-Kebir
 */

#include "agent.h"
#include "grgrid.h"
#include "serialization.h"

using namespace std;

const std::string Agent::_ClassName = "Agent";

unsigned long Agent::_nextID = 0;

// Needed for deserializing the model state.
// Avoids the calls to the random number generator in the normal constructor, allowing the random number generator
// to remain in synch after deserialization.
Agent::Agent()
    : _id(0)
	, _birthTime(-1)
	, _deathTime(-1)
	, _pos(-1, -1)
	, _trackMolecularDynamics(false)

	// TNFR components
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
    , _meanIL10R(0)

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

    // Valarrays for Numerical Methods
    , _initvector(0) // Used for numerical methods - current values of ODEs
    , _k1vector(0) // Forward Euler Solution (For both FE and RK)
    , _k2vector(0) // RK4 K2
    , _k3vector(0) // RK4 K3
    , _k4vector(0) // RK4 K4
    , _switchvector(0) // Valarray used to store solutions before reassigning to member variables

{
}

Agent::Agent(int birthtime, int deathtime, int row, int col

				//TNFR Components
				, Scalar meanTNFR1
				, Scalar stdTNFR1
				, Scalar meanTNFR2
				, Scalar stdTNFR2
				, Scalar kSynth
				, Scalar kTACE

				// IL10 components
				, Scalar iIL10R
				, Scalar stdIL10R
                , int odesize
			)
    : _id(0)
	, _birthTime(birthtime)
	, _deathTime(deathtime)
	, _pos(row, col)
	, _trackMolecularDynamics(false)

	// TNFR components
	, _mTNF(0.0)
	, _surfTNFR1(g_Rand.getReal(meanTNFR1 - stdTNFR1, meanTNFR1 + stdTNFR1))
	, _surfTNFR2(g_Rand.getReal(meanTNFR2 - stdTNFR2, meanTNFR2 + stdTNFR2))
	//	, _surfTNFR1(g_Rand.getLogNormal(meanTNFR1),_PARAM(stdTNFR1)))
	//	, _surfTNFR2(g_Rand.getLogNormal(meanTNFR2),_PARAM(stdTNFR2)))
	, _surfBoundTNFR1(0.0)
	, _surfBoundTNFR2(0.0)
	, _intBoundTNFR1(0.0)
	, _intBoundTNFR2(0.0)
	, _mTNFRNA(0.0)
	, _vTNFR1(_surfTNFR1 * _PARAM(PARAM_GR_K_T1))
	, _vTNFR2(_surfTNFR2 * _PARAM(PARAM_GR_K_T2))
	, _kSynth(kSynth)
	, _kTACE(kTACE)
	, _kmRNA(0.0)

	// IL10 components
	, _surfIL10R(g_Rand.getReal(iIL10R - stdIL10R, iIL10R + stdIL10R))
	, _vIL10R(_surfIL10R * _PARAM(PARAM_GR_I_K_T))
	, _surfBoundIL10R(0.0)
	, _kISynth(0.0)
    , _meanIL10R(iIL10R)

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

    // Valarrays for Numerical Methods
    , _initvector(0.0, odesize)
    , _k1vector(0.0, odesize)
    , _k2vector(0.0, odesize)
    , _k3vector(0.0, odesize)
    , _k4vector(0.0, odesize)
    , _switchvector(0.0, odesize)

{
	_id = createID();
}

Agent::~Agent()
{
}


void Agent::writeValarrayFromMembers(GrGrid& grid, valarray<double>& inputVector)
{
    // Defining the initial vector from the agent private variables
    // and getting soluble values from the grid
    int vectorsize = inputVector.size();
    
    if (vectorsize == 3) 
    {
        // sIL10
        inputVector[0] = grid.il10(_pos) / (NAV * VOL);
        // surfIL10R
        inputVector[1] = _surfIL10R;
        // surfBoundIL10R
        inputVector[2] = _surfBoundIL10R;
    }
    
    if (vectorsize > 3) 
    {
        
        // mTNFRNA
        inputVector[0] = _mTNFRNA;
        // mTNF
        inputVector[1] = _mTNF;
        // surfTNFR1
        inputVector[2] = _surfTNFR1;
        // surfTNFR2
        inputVector[3] = _surfTNFR2;
        // surfBoundTNFR1
        inputVector[4] = _surfBoundTNFR1;
        // surfBoundTNFR2
        inputVector[5] = _surfBoundTNFR2;
        // intBoundTNFR1 
        inputVector[6] = _intBoundTNFR1;
        // intBoundTNFR2
        inputVector[7] = _intBoundTNFR2;
        // sTNF
        inputVector[8] = grid.TNF(_pos) / (NAV * VOL);
        // shedTNFR2
        inputVector[9] = grid.shedTNFR2(_pos) / (NAV * VOL);
        
        if (vectorsize == 13) 
        {
            // sIL10
            inputVector[10] = grid.il10(_pos) / (NAV * VOL);
            // surfIL10R
            inputVector[11] = _surfIL10R;
            // surfBoundIL10R
            inputVector[12] = _surfBoundIL10R;
        }
        
        if (vectorsize == 38 || vectorsize == 41) 
        {
            inputVector[10] = _IKKKa;
            inputVector[11] = _IKKn; 
            inputVector[12] = _IKKa; 
            inputVector[13] = _IKKi;
            inputVector[14] = _IkBp; 
            inputVector[15] = _NFkB_IkBp; 
            inputVector[16] = _NFkBc; 
            inputVector[17] = _NFkBn; 
            inputVector[18] = _A20; 
            inputVector[19] = _A20t; 
            inputVector[20] = _IkB; 
            inputVector[21] = _IkBn; 
            inputVector[22] = _IkBt; 
            inputVector[23] = _NFkB_IkB; 
            inputVector[24] = _NFkB_IkBn; 
            inputVector[25] = _GA20; 
            inputVector[26] = _GIkB;  
            inputVector[27] = _GR; 
            inputVector[28] = _chemt; 
            inputVector[29] = _chem; 
            inputVector[30] = _TNFt; 
            inputVector[31] = _TNF; 
            inputVector[32] = _ACTt; 
            inputVector[33] = _ACT; 
            inputVector[34] = _normalizedACT; 
            inputVector[35] = _IAPt; 
            inputVector[36] = _IAP; 
            inputVector[37] = _normalizedIAP;
            
            if (vectorsize == 41) 
            {
                // sIL10
                inputVector[38] = grid.il10(_pos) / (NAV * VOL);
                // surfIL10R
                inputVector[39] = _surfIL10R;
                // surfBoundIL10R
                inputVector[40] = _surfBoundIL10R;
            }
        }
    }
}


void Agent::writeMembersFromValarray(GrGrid& grid, const valarray<double>& inputVector)
{
    // Defining the agent private variables and soluble variables from
    // the output vector of the numerical solver
    int vectorsize = inputVector.size();
    
    if (vectorsize == 3) 
    {
        // sIL10
        grid.setil10(_pos, (NAV * VOL * inputVector[0]));
        // surfIL10R
        _surfIL10R = inputVector[1];
        // surfBoundIL10R
        _surfBoundIL10R = inputVector[2];
    }
    
    if (vectorsize > 3) 
    {
        
        // mTNFRNA
        _mTNFRNA = inputVector[0];
        // mTNF
        _mTNF = inputVector[1];
        // surfTNFR1
        _surfTNFR1 = inputVector[2];
        // surfTNFR2
        _surfTNFR2 = inputVector[3]; 
        // surfBoundTNFR1
        _surfBoundTNFR1 = inputVector[4];
        // surfBoundTNFR2
        _surfBoundTNFR2 = inputVector[5]; 
        // intBoundTNFR1 
        _intBoundTNFR1 = inputVector[6]; 
        // intBoundTNFR2
        _intBoundTNFR2 = inputVector[7];
        // sTNF
        grid.setTNF(_pos, (NAV * VOL * inputVector[8]));
        // shedTNFR2
        grid.setshedTNFR2(_pos, (NAV * VOL * inputVector[9]));

        if (vectorsize == 13) 
        {
            // sIL10
            grid.setil10(_pos, (NAV * VOL * inputVector[10]));
            // surfIL10R
            _surfIL10R = inputVector[11];
            // surfBoundIL10R
            _surfBoundIL10R = inputVector[12];
        }
        
        if (vectorsize == 38 || vectorsize == 41) 
        {
            _IKKKa = inputVector[10];
            _IKKn = inputVector[11];
            _IKKa = inputVector[12];
            _IKKi = inputVector[13];
            _IkBp = inputVector[14];
            _NFkB_IkBp = inputVector[15];
            _NFkBc = inputVector[16];
            _NFkBn = inputVector[17];
            _A20 = inputVector[18];
            _A20t = inputVector[19];
            _IkB = inputVector[20];
            _IkBn = inputVector[21];
            _IkBt = inputVector[22];
            _NFkB_IkB = inputVector[23];
            _NFkB_IkBn = inputVector[24];
            _GA20 = inputVector[25];
            _GIkB = inputVector[26];
            _GR = inputVector[27];
            _chemt = inputVector[28];
            _chem = inputVector[29];
            _TNFt = inputVector[30];
            _TNF = inputVector[31];
            _ACTt = inputVector[32];
            _ACT = inputVector[33];
            _normalizedACT = inputVector[34];
            _IAPt = inputVector[35];
            _IAP = inputVector[36];
            _normalizedIAP = inputVector[37];
            
            grid.incCCL2(_pos, (_PARAM(PARAM_GR_e3Chem) * inputVector[29] * _PARAM(PARAM_GR_DT_MOLECULAR)));
            grid.incCCL5(_pos, (_PARAM(PARAM_GR_e3Chem) * inputVector[29] * _PARAM(PARAM_GR_DT_MOLECULAR)));
            grid.incCXCL9(_pos, (2 * _PARAM(PARAM_GR_e3Chem) * inputVector[29] * _PARAM(PARAM_GR_DT_MOLECULAR)));
            
            if (vectorsize == 41) 
            {
                // sIL10
                grid.setil10(_pos, (NAV * VOL * inputVector[38]));
                // surfIL10R
                _surfIL10R = inputVector[39];
                // surfBoundIL10R
                _surfBoundIL10R = inputVector[40];
            }
        }
    }
}







void Agent::derivativeTNFandIL10(const valarray<double>& vecread, valarray<double>& vecwrite, double dt, GrGrid& grid)
{
    
}

void Agent::derivativeTNFandNFKB(const valarray<double>& vecread, valarray<double>& vecwrite, double dt, GrGrid& grid)
{
    
}

void Agent::derivativeTNFandIL10andNFKB(const valarray<double>& vecread, valarray<double>& vecwrite, double dt, GrGrid& grid)
{
    
}

void Agent::derivativeTNF(const valarray<double>& vecread, valarray<double>& vecwrite, double dt, GrGrid& grid)
{
    assert(vecread.size()==10);
    
    double koff1 = _PARAM(PARAM_GR_K_ON1) * _PARAM(PARAM_GR_KD1);
	double koff2 = _PARAM(PARAM_GR_K_ON2) * _PARAM(PARAM_GR_KD2);
    double il10 = grid.il10(_pos) /(NAV * VOL);
    double IkmRNA;
    
    // modulate TNF parameters based on equilibrium IL10 receptor equations since the il10 molecular equations are not active
    double eqsurfBoundIL10R = (il10 * _meanIL10R) / (_PARAM(PARAM_GR_I_KD) + il10);
    
    if (_kmRNA > 0) {
        IkmRNA = _kmRNA * ((_kSynth/_kmRNA)+ ((1.0 - (_kSynth/_kmRNA))/(1.0 + pow(2.7183, ((eqsurfBoundIL10R - _PARAM(PARAM_GR_LINK_RNA_GAMMA))/_PARAM(PARAM_GR_LINK_RNA_DELTA))))));
    }
    else
    {
        IkmRNA = 0.0;
    }
    
    // TNF Ordinary Differential Equations
    // mTNFRNA
    vecwrite[0] = (IkmRNA - _PARAM(PARAM_GR_K_TRANS) * vecread[0]) * dt;
    // mTNF
    vecwrite[1] = (_PARAM(PARAM_GR_K_TRANS) * vecread[0] - _kTACE * vecread[1]) * dt;
    // surfTNFR1
    vecwrite[2] = (_vTNFR1 - _PARAM(PARAM_GR_K_ON1) * vecread[8] * vecread[2] + koff1 * vecread[4] - _PARAM(PARAM_GR_K_T1) * vecread[2] + _PARAM(PARAM_GR_K_REC1) * vecread[6]) * dt;
    // surfTNFR2
    vecwrite[3] = (_vTNFR2 - _PARAM(PARAM_GR_K_ON2) * vecread[8] * vecread[3] + koff2 * vecread[5] - _PARAM(PARAM_GR_K_T2) * vecread[3] + _PARAM(PARAM_GR_K_REC2) * vecread[7]) * dt;
    // surfBoundTNFR1
    vecwrite[4] = (_PARAM(PARAM_GR_K_ON1) * vecread[8] * vecread[2] - koff1 * vecread[4] - _PARAM(PARAM_GR_K_INT1) * vecread[4]) * dt;
    // surfBoundTNFR2
    vecwrite[5] = (_PARAM(PARAM_GR_K_ON2) * vecread[8] * vecread[3] - koff2 * vecread[5] - _PARAM(PARAM_GR_K_INT2) * vecread[5] - _PARAM(PARAM_GR_K_SHED) * vecread[5]) * dt;
    // intBoundTNFR1 
    vecwrite[6] = (_PARAM(PARAM_GR_K_INT1) * vecread[4] - _PARAM(PARAM_GR_K_DEG1) * vecread[6] - _PARAM(PARAM_GR_K_REC1) * vecread[6]) * dt;
    // intBoundTNFR2
    vecwrite[7] = (_PARAM(PARAM_GR_K_INT2) * vecread[5] - _PARAM(PARAM_GR_K_DEG2) * vecread[7] - _PARAM(PARAM_GR_K_REC2) * vecread[7]) * dt;
    // sTNF
    vecwrite[8] = ((DENSITY/NAV) * (_kTACE * vecread[1] - _PARAM(PARAM_GR_K_ON1) * vecread[8] * vecread[2] + koff1 * vecread[4] - _PARAM(PARAM_GR_K_ON2) * vecread[8] * vecread[3] + koff2 * vecread[5])) * dt;
    // shedTNFR2
    vecwrite[9] = ((DENSITY/NAV) * _PARAM(PARAM_GR_K_SHED) * vecread[5]) * dt;
    
}

void Agent::solveTNF(GrGrid& grid, double dt)
{
	double koff1 = _PARAM(PARAM_GR_K_ON1) * _PARAM(PARAM_GR_KD1);
	double koff2 = _PARAM(PARAM_GR_K_ON2) * _PARAM(PARAM_GR_KD2);

	double tnf = grid.TNF(_pos) / (NAV * VOL);
	double shedtnfr2 = grid.shedTNFR2(_pos) / (NAV * VOL);
    double il10 = grid.il10(_pos) /(NAV * VOL);

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
    eqsurfBoundIL10R = (il10 * _meanIL10R) / (_PARAM(PARAM_GR_I_KD) + il10);

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
	dsTNF = ((DENSITY/NAV) * (_kTACE * _mTNF - _PARAM(PARAM_GR_K_ON1) * tnf * _surfTNFR1 + koff1 * _surfBoundTNFR1 - _PARAM(PARAM_GR_K_ON2) * tnf * _surfTNFR2 + koff2 * _surfBoundTNFR2)) * dt;
	dshedTNFR2 = ((DENSITY/NAV) * _PARAM(PARAM_GR_K_SHED) * _surfBoundTNFR2) * dt;

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

	grid.setTNF(_pos, (NAV * VOL * tnf));
	grid.setshedTNFR2(_pos, (NAV * VOL * shedtnfr2));
	if (_mTNF < 0 || _surfTNFR1 < 0 || _surfBoundTNFR1 < 0 || _surfTNFR2 < 0 || _surfBoundTNFR2 < 0 || _mTNFRNA < 0)
		std::cout << "Error: Negative Value of Species in TNF/TNFR dynamics" << std::endl;

    //cout << "Debug: Running TNF dynamics" << std::endl;
}

void Agent::derivativeIL10(const valarray<double>& vecread, valarray<double>& vecwrite, double dt, GrGrid& grid)
{
    // IL10 Ordinary Differential Equations
    // sIL10
    vecwrite[0] = (((DENSITY/NAV) * _kISynth) + ((DENSITY/NAV) * (_PARAM(PARAM_GR_I_K_OFF) * vecread[2] - _PARAM(PARAM_GR_I_K_ON) * vecread[1] * vecread[0]))) * dt;
    // surfIL10R
    vecwrite[1] = (_vIL10R - _PARAM(PARAM_GR_I_K_ON) * vecread[1] * vecread[0] + _PARAM(PARAM_GR_I_K_OFF) * vecread[2] - _PARAM(PARAM_GR_I_K_T) * vecread[1]) * dt;
    // surfBoundIL10R
    vecwrite[2] = (_PARAM(PARAM_GR_I_K_ON) * vecread[1] * vecread[0] - _PARAM(PARAM_GR_I_K_OFF) * vecread[2] - _PARAM(PARAM_GR_I_K_INT) * vecread[2]) * dt;
}

void Agent::solveIL10(GrGrid& grid, double dt)
{
    double il10 = grid.il10(_pos) / (NAV * VOL);

    double dsIL10;
	double dsurfIL10R;
	double dsurfBoundIL10R;

    // IL10 differential equations
    dsIL10 = (DENSITY/NAV) * _kISynth + ((DENSITY/NAV) * (_PARAM(PARAM_GR_I_K_OFF) * _surfBoundIL10R - _PARAM(PARAM_GR_I_K_ON) * _surfIL10R * il10)) * dt;
	dsurfIL10R = (_vIL10R - _PARAM(PARAM_GR_I_K_ON) * _surfIL10R * il10 + _PARAM(PARAM_GR_I_K_OFF) * _surfBoundIL10R - _PARAM(PARAM_GR_I_K_T) * _surfIL10R) * dt;
	dsurfBoundIL10R = (_PARAM(PARAM_GR_I_K_ON) * _surfIL10R * il10 - _PARAM(PARAM_GR_I_K_OFF) * _surfBoundIL10R - _PARAM(PARAM_GR_I_K_INT) * _surfBoundIL10R) * dt;
    // end of IL10 differential equations

    // update il10 variables
    _surfIL10R += dsurfIL10R;
    _surfBoundIL10R += dsurfBoundIL10R;
    il10 += dsIL10;

    grid.setil10(_pos, (NAV * VOL * il10));

    if (_surfIL10R < 0 || _surfBoundIL10R < 0)
        std::cout << "Error: Negative value of species in IL10/IL10R dynamics" << std::endl;

     //cout << "Debug: Running IL10 dynamics" << std::endl;

}

void Agent::solveTNFandIL10(GrGrid& grid, double dt)
{
	double koff1 = _PARAM(PARAM_GR_K_ON1) * _PARAM(PARAM_GR_KD1);
	double koff2 = _PARAM(PARAM_GR_K_ON2) * _PARAM(PARAM_GR_KD2);
	
	double tnf = grid.TNF(_pos) / (NAV * VOL);
	double shedtnfr2 = grid.shedTNFR2(_pos) / (NAV * VOL);
    double il10 = grid.il10(_pos) / (NAV * VOL);
	
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
    
    if (_kmRNA > 0)
    {
        IkmRNA = _kmRNA * ((_kSynth/_kmRNA) + ((1.0 - (_kSynth/_kmRNA))/(1.0 + pow(2.7183, ((_surfBoundIL10R - _PARAM(PARAM_GR_LINK_RNA_GAMMA))/_PARAM(PARAM_GR_LINK_RNA_DELTA))))));
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
	dsTNF = ((DENSITY/NAV) * (_kTACE * _mTNF - _PARAM(PARAM_GR_K_ON1) * tnf * _surfTNFR1 + koff1 * _surfBoundTNFR1 - _PARAM(PARAM_GR_K_ON2) * tnf * _surfTNFR2 + koff2 * _surfBoundTNFR2)) * dt;
	dshedTNFR2 = ((DENSITY/NAV) * _PARAM(PARAM_GR_K_SHED) * _surfBoundTNFR2) * dt;
	// end of TNF differential equations
    
    // IL10 differential equations
    dsIL10 = (DENSITY/NAV) * _kISynth + ((DENSITY/NAV) * (_PARAM(PARAM_GR_I_K_OFF) * _surfBoundIL10R - _PARAM(PARAM_GR_I_K_ON) * _surfIL10R * il10)) * dt;
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

	grid.setTNF(_pos, (NAV * VOL * tnf));
	grid.setshedTNFR2(_pos, (NAV * VOL * shedtnfr2));
	grid.setil10(_pos, (NAV * VOL * il10));
	
    
	if (_mTNF < 0 || _surfTNFR1 < 0 || _surfBoundTNFR1 < 0 || _surfTNFR2 < 0 || _surfBoundTNFR2 < 0 || _mTNFRNA < 0)
		std::cout << "Error: Negative Value of Species in TNF/TNFR dynamics" << std::endl;
    
    if (_surfIL10R < 0 || _surfBoundIL10R < 0)
        std::cout << "Error: Negative value of species in IL10/IL10R dynamics" << std::endl;
    
}

void Agent::solveNFkBandTNF(GrGrid& grid, double dt)
{
	double koff1 = _PARAM(PARAM_GR_K_ON1) * _PARAM(PARAM_GR_KD1);
	double koff2 = _PARAM(PARAM_GR_K_ON2) * _PARAM(PARAM_GR_KD2);
	
	double tnf = grid.TNF(_pos) / (NAV * VOL);
	double shedtnfr2 = grid.shedTNFR2(_pos) / (NAV * VOL);
	double il10 = grid.il10(_pos) /(NAV * VOL);
    
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
	dsurfTNFR1 = (_vTNFR1 - _PARAM(PARAM_GR_K_ON1) * tnf * _surfTNFR1 + (koff1+KDEG) * _surfBoundTNFR1 - _PARAM(PARAM_GR_K_T1) * _surfTNFR1 + _PARAM(PARAM_GR_K_REC1) * _intBoundTNFR1) * dt;
	dsurfTNFR2 = (_vTNFR2 - _PARAM(PARAM_GR_K_ON2) * tnf * _surfTNFR2 + (koff2+KDEG) * _surfBoundTNFR2 - _PARAM(PARAM_GR_K_T2) * _surfTNFR2 + _PARAM(PARAM_GR_K_REC2) * _intBoundTNFR2) * dt;
	dsurfBoundTNFR1 = (_PARAM(PARAM_GR_K_ON1) * tnf * _surfTNFR1 - (koff1+KDEG) * _surfBoundTNFR1 - _PARAM(PARAM_GR_K_INT1) * _surfBoundTNFR1) * dt;
	dsurfBoundTNFR2 = (_PARAM(PARAM_GR_K_ON2) * tnf * _surfTNFR2 - (koff2+KDEG) * _surfBoundTNFR2 - _PARAM(PARAM_GR_K_INT2) * _surfBoundTNFR2 - _PARAM(PARAM_GR_K_SHED) * _surfBoundTNFR2) * dt;
	dintBoundTNFR1 = (_PARAM(PARAM_GR_K_INT1) * _surfBoundTNFR1 - _PARAM(PARAM_GR_K_DEG1) * _intBoundTNFR1 - _PARAM(PARAM_GR_K_REC1) * _intBoundTNFR1) * dt;
	dintBoundTNFR2 = (_PARAM(PARAM_GR_K_INT2) * _surfBoundTNFR2 - _PARAM(PARAM_GR_K_DEG2) * _intBoundTNFR2 - _PARAM(PARAM_GR_K_REC2) * _intBoundTNFR2) * dt;
	dsTNF = ((DENSITY/NAV) * (_kTACE * _mTNF - _PARAM(PARAM_GR_K_ON1) * tnf * _surfTNFR1 + koff1 * _surfBoundTNFR1 - _PARAM(PARAM_GR_K_ON2) * tnf * _surfTNFR2 + koff2 * _surfBoundTNFR2)) * dt;
	dshedTNFR2 = ((DENSITY/NAV) * _PARAM(PARAM_GR_K_SHED) * _surfBoundTNFR2) * dt;
	
	// NF-kB dynamics model equations
	dIKKKa = (_PARAM(PARAM_GR_ka)*_surfBoundTNFR1*(_PARAM(PARAM_GR_KN)-_IKKKa)*_PARAM(PARAM_GR_kA20)/(_PARAM(PARAM_GR_kA20)+_A20)-_PARAM(PARAM_GR_ki)*_IKKKa) * dt;
	dIKKn = (_PARAM(PARAM_GR_k4)*(_PARAM(PARAM_GR_KNN)-_IKKn-_IKKa-_IKKi)-_PARAM(PARAM_GR_k1)*pow(_IKKKa,2.0)*_IKKn) * dt;
	dIKKa = (_PARAM(PARAM_GR_k1)*pow(_IKKKa,2.0)*_IKKn-_PARAM(PARAM_GR_k3)*_IKKa*(_PARAM(PARAM_GR_k2)+_A20)/_PARAM(PARAM_GR_k2)) * dt;
	dIKKi = (_PARAM(PARAM_GR_k3)*_IKKa*(_PARAM(PARAM_GR_k2)+_A20)/_PARAM(PARAM_GR_k2)-_PARAM(PARAM_GR_k4)*_IKKi) * dt;
	dIkBp = (_PARAM(PARAM_GR_a2)*_IKKa*_IkB-_PARAM(PARAM_GR_tp)*_IkBp) * dt;
	dNFkB_IkBp = (_PARAM(PARAM_GR_a3)*_IKKa*_NFkB_IkB-_PARAM(PARAM_GR_tp)*_NFkB_IkBp) * dt;
	dNFkBc = (_PARAM(PARAM_GR_c6a)*_NFkB_IkB-_PARAM(PARAM_GR_a1)*_NFkBc*_IkB+_PARAM(PARAM_GR_tp)*_NFkB_IkBp-_PARAM(PARAM_GR_i1)*_NFkBc) * dt;
	dNFkBn = (_PARAM(PARAM_GR_i1)*_NFkBc-_PARAM(PARAM_GR_a1)*KV*_IkBn*_NFkBn) * dt;
	dA20 = (_PARAM(PARAM_GR_c4)*_A20t-_PARAM(PARAM_GR_c5)*_A20) * dt;
	dA20t = (_PARAM(PARAM_GR_c1)*_GA20-_PARAM(PARAM_GR_c3)*_A20t) * dt;
	dIkB = (-_PARAM(PARAM_GR_a2)*_IKKa*_IkB-_PARAM(PARAM_GR_a1)*_IkB*_NFkBc+_PARAM(PARAM_GR_c4)*_IkBt-_PARAM(PARAM_GR_c5a)*_IkB-_PARAM(PARAM_GR_i1a)*_IkB+_PARAM(PARAM_GR_e1a)*_IkBn) * dt;
	dIkBn = (-_PARAM(PARAM_GR_a1)*KV*_IkBn*_NFkBn+_PARAM(PARAM_GR_i1a)*_IkB-_PARAM(PARAM_GR_e1a)*_IkBn) * dt;
	dIkBt = (_PARAM(PARAM_GR_c1)*_GIkB-_PARAM(PARAM_GR_c3)*_IkBt) * dt;
	dNFkB_IkB = (_PARAM(PARAM_GR_a1)*_IkB*_NFkBc-_PARAM(PARAM_GR_c6a)*_NFkB_IkB-_PARAM(PARAM_GR_a3)*_IKKa*_NFkB_IkB+_PARAM(PARAM_GR_e2a)*_NFkB_IkBn) * dt;
	dNFkB_IkBn = (_PARAM(PARAM_GR_a1)*KV*_IkBn*_NFkBn-_PARAM(PARAM_GR_e2a)*_NFkB_IkBn) * dt;
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
	
	grid.setTNF(_pos, (NAV * VOL * tnf));
	grid.setshedTNFR2(_pos, (NAV * VOL * shedtnfr2));
	
	// secrete chemokines
	grid.incCCL2(_pos, (_PARAM(PARAM_GR_e3Chem) * _chem * dt));
	grid.incCCL5(_pos, (_PARAM(PARAM_GR_e3Chem) * _chem * dt));
	grid.incCXCL9(_pos, (2 * _PARAM(PARAM_GR_e3Chem) * _chem * dt));
	
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

void Agent::solveTNFandIL10andNFkB(GrGrid& grid, double dt)
{
	double koff1 = _PARAM(PARAM_GR_K_ON1) * _PARAM(PARAM_GR_KD1);
	double koff2 = _PARAM(PARAM_GR_K_ON2) * _PARAM(PARAM_GR_KD2);
	
	double tnf = grid.TNF(_pos) / (NAV * VOL);
	double shedtnfr2 = grid.shedTNFR2(_pos) / (NAV * VOL);
    double il10 = grid.il10(_pos) / (NAV * VOL);
	
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
	dsurfTNFR1 = (_vTNFR1 - _PARAM(PARAM_GR_K_ON1) * tnf * _surfTNFR1 + (koff1+KDEG) * _surfBoundTNFR1 - _PARAM(PARAM_GR_K_T1) * _surfTNFR1 + _PARAM(PARAM_GR_K_REC1) * _intBoundTNFR1) * dt;
	dsurfTNFR2 = (_vTNFR2 - _PARAM(PARAM_GR_K_ON2) * tnf * _surfTNFR2 + (koff2+KDEG) * _surfBoundTNFR2 - _PARAM(PARAM_GR_K_T2) * _surfTNFR2 + _PARAM(PARAM_GR_K_REC2) * _intBoundTNFR2) * dt;
	dsurfBoundTNFR1 = (_PARAM(PARAM_GR_K_ON1) * tnf * _surfTNFR1 - (koff1+KDEG) * _surfBoundTNFR1 - _PARAM(PARAM_GR_K_INT1) * _surfBoundTNFR1) * dt;
	dsurfBoundTNFR2 = (_PARAM(PARAM_GR_K_ON2) * tnf * _surfTNFR2 - (koff2+KDEG) * _surfBoundTNFR2 - _PARAM(PARAM_GR_K_INT2) * _surfBoundTNFR2 - _PARAM(PARAM_GR_K_SHED) * _surfBoundTNFR2) * dt;
	dintBoundTNFR1 = (_PARAM(PARAM_GR_K_INT1) * _surfBoundTNFR1 - _PARAM(PARAM_GR_K_DEG1) * _intBoundTNFR1 - _PARAM(PARAM_GR_K_REC1) * _intBoundTNFR1) * dt;
	dintBoundTNFR2 = (_PARAM(PARAM_GR_K_INT2) * _surfBoundTNFR2 - _PARAM(PARAM_GR_K_DEG2) * _intBoundTNFR2 - _PARAM(PARAM_GR_K_REC2) * _intBoundTNFR2) * dt;
	dsTNF = ((DENSITY/NAV) * (_kTACE * _mTNF - _PARAM(PARAM_GR_K_ON1) * tnf * _surfTNFR1 + koff1 * _surfBoundTNFR1 - _PARAM(PARAM_GR_K_ON2) * tnf * _surfTNFR2 + koff2 * _surfBoundTNFR2)) * dt;
	dshedTNFR2 = ((DENSITY/NAV) * _PARAM(PARAM_GR_K_SHED) * _surfBoundTNFR2) * dt;
	
	// NF-kB dynamics model equations
	dIKKKa = (_PARAM(PARAM_GR_ka)*_surfBoundTNFR1*(_PARAM(PARAM_GR_KN)-_IKKKa)*_PARAM(PARAM_GR_kA20)/(_PARAM(PARAM_GR_kA20)+_A20)-_PARAM(PARAM_GR_ki)*_IKKKa) * dt;
	dIKKn = (_PARAM(PARAM_GR_k4)*(_PARAM(PARAM_GR_KNN)-_IKKn-_IKKa-_IKKi)-_PARAM(PARAM_GR_k1)*pow(_IKKKa,2.0)*_IKKn) * dt;
	dIKKa = (_PARAM(PARAM_GR_k1)*pow(_IKKKa,2.0)*_IKKn-_PARAM(PARAM_GR_k3)*_IKKa*(_PARAM(PARAM_GR_k2)+_A20)/_PARAM(PARAM_GR_k2)) * dt;
	dIKKi = (_PARAM(PARAM_GR_k3)*_IKKa*(_PARAM(PARAM_GR_k2)+_A20)/_PARAM(PARAM_GR_k2)-_PARAM(PARAM_GR_k4)*_IKKi) * dt;
	dIkBp = (_PARAM(PARAM_GR_a2)*_IKKa*_IkB-_PARAM(PARAM_GR_tp)*_IkBp) * dt;
	dNFkB_IkBp = (_PARAM(PARAM_GR_a3)*_IKKa*_NFkB_IkB-_PARAM(PARAM_GR_tp)*_NFkB_IkBp) * dt;
	dNFkBc = (_PARAM(PARAM_GR_c6a)*_NFkB_IkB-_PARAM(PARAM_GR_a1)*_NFkBc*_IkB+_PARAM(PARAM_GR_tp)*_NFkB_IkBp-_PARAM(PARAM_GR_i1)*_NFkBc) * dt;
	dNFkBn = (_PARAM(PARAM_GR_i1)*_NFkBc-_PARAM(PARAM_GR_a1)*KV*_IkBn*_NFkBn) * dt;
	dA20 = (_PARAM(PARAM_GR_c4)*_A20t-_PARAM(PARAM_GR_c5)*_A20) * dt;
	dA20t = (_PARAM(PARAM_GR_c1)*_GA20-_PARAM(PARAM_GR_c3)*_A20t) * dt;
	dIkB = (-_PARAM(PARAM_GR_a2)*_IKKa*_IkB-_PARAM(PARAM_GR_a1)*_IkB*_NFkBc+_PARAM(PARAM_GR_c4)*_IkBt-_PARAM(PARAM_GR_c5a)*_IkB-_PARAM(PARAM_GR_i1a)*_IkB+_PARAM(PARAM_GR_e1a)*_IkBn) * dt;
	dIkBn = (-_PARAM(PARAM_GR_a1)*KV*_IkBn*_NFkBn+_PARAM(PARAM_GR_i1a)*_IkB-_PARAM(PARAM_GR_e1a)*_IkBn) * dt;
	dIkBt = (_PARAM(PARAM_GR_c1)*_GIkB-_PARAM(PARAM_GR_c3)*_IkBt) * dt;
	dNFkB_IkB = (_PARAM(PARAM_GR_a1)*_IkB*_NFkBc-_PARAM(PARAM_GR_c6a)*_NFkB_IkB-_PARAM(PARAM_GR_a3)*_IKKa*_NFkB_IkB+_PARAM(PARAM_GR_e2a)*_NFkB_IkBn) * dt;
	dNFkB_IkBn = (_PARAM(PARAM_GR_a1)*KV*_IkBn*_NFkBn-_PARAM(PARAM_GR_e2a)*_NFkB_IkBn) * dt;
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
    dsIL10 = (DENSITY/NAV) * _kISynth + ((DENSITY/NAV) * (_PARAM(PARAM_GR_I_K_OFF) * _surfBoundIL10R - _PARAM(PARAM_GR_I_K_ON) * _surfIL10R * il10)) * dt;
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
    
	grid.setTNF(_pos, (NAV * VOL * tnf));
	grid.setshedTNFR2(_pos, (NAV * VOL * shedtnfr2));
	grid.setil10(_pos, (NAV * VOL * il10));
    
	// secrete chemokines
	grid.incCCL2(_pos, (_PARAM(PARAM_GR_e3Chem) * _chem * dt));
	grid.incCCL5(_pos, (_PARAM(PARAM_GR_e3Chem) * _chem * dt));
	grid.incCXCL9(_pos,  (2 * _PARAM(PARAM_GR_e3Chem) * _chem * dt));
	
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

void Agent::solveNFkBODEsEquilibrium(double dt)
{
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
	dNFkBn = (_PARAM(PARAM_GR_i1)*_NFkBc-_PARAM(PARAM_GR_a1)*KV*_IkBn*_NFkBn) * dt;
	dA20 = (_PARAM(PARAM_GR_c4)*_A20t-_PARAM(PARAM_GR_c5)*_A20) * dt;
	dA20t = (_PARAM(PARAM_GR_c1)*_GA20-_PARAM(PARAM_GR_c3)*_A20t) * dt;
	dIkB = (-_PARAM(PARAM_GR_a2)*_IKKa*_IkB-_PARAM(PARAM_GR_a1)*_IkB*_NFkBc+_PARAM(PARAM_GR_c4)*_IkBt-_PARAM(PARAM_GR_c5a)*_IkB-_PARAM(PARAM_GR_i1a)*_IkB+_PARAM(PARAM_GR_e1a)*_IkBn) * dt;
	dIkBn = (-_PARAM(PARAM_GR_a1)*KV*_IkBn*_NFkBn+_PARAM(PARAM_GR_i1a)*_IkB-_PARAM(PARAM_GR_e1a)*_IkBn) * dt;
	dIkBt = (_PARAM(PARAM_GR_c1)*_GIkB-_PARAM(PARAM_GR_c3)*_IkBt) * dt;
	dNFkB_IkB = (_PARAM(PARAM_GR_a1)*_IkB*_NFkBc-_PARAM(PARAM_GR_c6a)*_NFkB_IkB-_PARAM(PARAM_GR_a3)*_IKKa*_NFkB_IkB+_PARAM(PARAM_GR_e2a)*_NFkB_IkBn) * dt;
	dNFkB_IkBn = (_PARAM(PARAM_GR_a1)*KV*_IkBn*_NFkBn-_PARAM(PARAM_GR_e2a)*_NFkB_IkBn) * dt;
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

void Agent::solveDegradation(GrGrid& grid, double dt, bool tnfrDynamics, bool il10rDynamics, Scalar meanTNFR1, Scalar iIL10R)
{
    if (!tnfrDynamics) {

        // simulate the effect of TNF internalization by cells in the form of degradation. Only for TNF
        double dtnf;
        double tnf = grid.TNF(_pos);
        dtnf = -_PARAM(PARAM_GR_K_INT1) * (tnf / (tnf + _PARAM(PARAM_GR_KD1) * NAV * VOL)) * meanTNFR1 * dt * 0.4;
        grid.incTNF(_pos, dtnf);
    }

    if (!il10rDynamics) {

        double dil10;
        double il10 = grid.il10(_pos);

        // simulate the effect of IL10 internalization in the form of degradation. Only for IL10
        dil10 = -_PARAM(PARAM_GR_I_K_INT) * (il10 / (il10 + _PARAM(PARAM_GR_I_KD) * NAV * VOL)) * iIL10R * dt * _PARAM(PARAM_GR_I_MOD);
        grid.incil10(_pos, dil10);

    }

}

Pos Agent::moveAgent(GrGrid& grid, bool ccl2, bool ccl5, bool cxcl9, bool attractant, double bonusFactor)
{
	int DestinationOrdinal = getDestinationOrdinal(grid, ccl2, ccl5, cxcl9, attractant, bonusFactor);
	return compartmentOrdinalToCoordinates(DestinationOrdinal, grid.getRange());
}

int Agent::getDestinationOrdinal(GrGrid& grid, bool ccl2, bool ccl5, bool cxcl9, bool attractant, double bonusFactor)
{

	double min = _PARAM(PARAM_GR_MIN_CHEMOTAXIS);
	double max = _PARAM(PARAM_GR_MAX_CHEMOTAXIS);

	bool ccl2Switch =  min < grid.CCL2(_pos) && grid.CCL2(_pos) < max;
	bool ccl5Switch =  min < grid.CCL5(_pos) && grid.CCL5(_pos) < max;
	bool cxcl9Switch = min < grid.CXCL9(_pos)  && grid.CXCL9(_pos) < max;

	Scalar prob[MOORE_COUNT] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

  int _row(GETROW(_pos)), _col(GETCOL(_pos));
	int k = 0;
  for (int i = -1; i <= 1; i++)
  {
    for (int j = -1; j <= 1; j++)
    {
      if (ccl2 && ccl2Switch)
        prob[k] += grid.CCL2(Pos(grid.mod_row(_row + i), grid.mod_col(_col + j)));
      if (ccl5 && ccl5Switch)
        prob[k] += grid.CCL5(Pos(grid.mod_row(_row + i), grid.mod_col(_col + j)));
      if (cxcl9 && cxcl9Switch)
        prob[k] += grid.CXCL9(Pos(grid.mod_row(_row + i), grid.mod_col(_col + j)));
      if (attractant)
        prob[k] += grid.macAttractant(Pos(grid.mod_row(_row + i), grid.mod_col(_col + j)));
      k++;
    }
  }
	
	// multiply highest probability with bonusFactor,
	// and while we are at it, determine the sum
	k = 0;
	double sum = 0;
	for (int i = 0; i < MOORE_COUNT; i++)
	{
		if (prob[i] > prob[k])
			k = i;
		sum += prob[i];
	}
	sum += (bonusFactor - 1) * prob[k];
	prob[k] *= bonusFactor;

	if (sum > 0.0)
	{
		// normalize
		for (int i = 0; i < MOORE_COUNT; i++)
			prob[i] /= sum;

		// compute cumulative array
		double cumProb[MOORE_COUNT];
    cumProb[0] = prob[0];
		for (int i = 1; i < MOORE_COUNT; i++)
		{
			cumProb[i] = (cumProb[i - 1] + prob[i]);
		}

		// linear search
		double r = g_Rand.getReal();
		for (k = 0; k < 9 && cumProb[k] < r; k++);
	}
	else
	{
		// prob[i] = 0 for all i.
		// Pick from the neighbors with equal probability.
		k = g_Rand.getInt(MOORE_COUNT, 0);
	}

	/**
	 * 0 1 2
	 * 3 4 5
	 * 6 7 8
	 */
  assert(0 <= k && k < MOORE_COUNT);
	return k;
}

Pos Agent::compartmentOrdinalToCoordinates(int ordinal, const Pos& dim) const
{

	/**
	 * The possible values for the ordinal are:
	 *
	 *  0  1  2
	 *  3  4  5
	 *  6  7  8
	 */
	int dRow = ((ordinal / 3) % 3) - 1;
	int dCol = ordinal % 3 - 1;
	int newRow = (_pos.x + dRow + dim.x) % dim.x;
	int newCol = (_pos.y + dCol + dim.y) % dim.y;

	return Pos(newRow, newCol);

}

void Agent::classSerialize(std::ostream& out)
{
	assert(out.good());
	Serialization::writeHeader(out, Agent::_ClassName);
	out << _nextID << std::endl;
	Serialization::writeFooter(out, Agent::_ClassName);
}

void Agent::classDeserialize(std::istream& in)
{
	assert(in.good());

	if (!Serialization::readHeader(in, Agent::_ClassName))
	{
		exit(1);
	}

	in >> _nextID;

	if (!Serialization::readFooter(in, Agent::_ClassName))
	{
		exit(1);
	}
}

void Agent::serialize(std::ostream& out) const
{
	assert(out.good());

	Serialization::writeHeader(out, Agent::_ClassName);

	out << _id << std::endl;
	out << _birthTime << std::endl;
	out << _deathTime << std::endl;
	out << _pos.x << std::endl;
	out << _pos.y << std::endl;
	out << _trackMolecularDynamics << std::endl;

	// TNF associated attributes
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

    // IL10 associated attributes
    out << _surfIL10R << std::endl;
    out << _vIL10R << std::endl;
    out << _surfBoundIL10R << std::endl;
    out << _kISynth << std::endl;
    out << _meanIL10R <<std::endl;
    

    // NF-kB signaling pathway components
	out << _IKKKa << std::endl; // (IKKK in active state)
	out << _IKKn << std::endl; // (IKK in neutral state)
	out << _IKKa << std::endl; // (IKK in the active state)
	out << _IKKi << std::endl; // (IKK in inactive state)
	out << _IkBp << std::endl; // (Phospho-IkB)
	out << _NFkB_IkBp << std::endl; // %NFkB|IkBp
	out << _NFkBc << std::endl; // cytoplasmic NFkB
	out << _NFkBn << std::endl; // nucluar NFkB
	out << _A20 << std::endl;
	out << _A20t << std::endl; // A20 transcript
	out << _IkB << std::endl;
	out << _IkBn << std::endl; // nucluar IkB
	out << _IkBt << std::endl; // IkB trancript
	out << _NFkB_IkB << std::endl;
	out << _NFkB_IkBn << std::endl;
	out << _GA20 << std::endl; // (The state of A20)
	out << _GIkB << std::endl; // (The state of IkB)
	out << _GR << std::endl; // (The state of reporter genes)
	out << _c1rrChemTNF << std::endl; // NF-kB independent rate of TNF/chemokine mRNA synthesis  
	out << _c1rChem << std::endl;
	out << _c1rTNF << std::endl;
	out << _chemt << std::endl; // generic chemokine transcript
	out << _chem << std::endl; // intracellular generic chemokine protein
	out << _TNFt << std::endl; // TNF transcript
	out << _TNF << std::endl; // intracellular TNF 
	out << _ACTt << std::endl; // transcript of macrophage activating molecules
	out << _ACT << std::endl; // mac-activation molecules
	out << _normalizedACT << std::endl;
	out << _IAPt << std::endl; // transcript of IAP (inhibitor of apoptosis)
	out << _IAP << std::endl; 
	out << _normalizedIAP << std::endl;

    out << _initvector.size() << endl;
    
    for (size_t jj = 0; jj < _initvector.size(); jj++) 
    {
        out << _initvector[jj] << std::endl;
        out << _k1vector[jj] << std::endl;
        out << _k2vector[jj] << std::endl;
        out << _k3vector[jj] << std::endl;
        out << _k4vector[jj] << std::endl;
        out << _switchvector[jj] << std::endl;
    }

	Serialization::writeFooter(out, Agent::_ClassName);
}

void Agent::deserialize(std::istream& in)
{
	assert(in.good());

	if (!Serialization::readHeader(in, Agent::_ClassName))
	{
		exit(1);
	}

	in >> _id;
	in >> _birthTime;
	in >> _deathTime;
	in >> _pos.x;
	in >> _pos.y;
	in >> _trackMolecularDynamics;

	// TNF associated attributes
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

    // IL10 associated attributes
    in >> _surfIL10R;
    in >> _vIL10R;
    in >> _surfBoundIL10R;
    in >> _kISynth;
    in >> _meanIL10R;


    // NF-kB signaling pathway components
	in >> _IKKKa; // (IKKK in active state)
	in >> _IKKn; // (IKK in neutral state)
	in >> _IKKa; // (IKK in the active state)
	in >> _IKKi; // (IKK in inactive state)
	in >> _IkBp; // (Phospho-IkB)
	in >> _NFkB_IkBp; // %NFkB|IkBp
	in >> _NFkBc; // cytoplasmic NFkB
	in >> _NFkBn; // nucluar NFkB
	in >> _A20;
	in >> _A20t; // A20 transcript
	in >> _IkB;
	in >> _IkBn; // nucluar IkB
	in >> _IkBt; // IkB trancript
	in >> _NFkB_IkB;
	in >> _NFkB_IkBn;
	in >> _GA20; // (The state of A20)
	in >> _GIkB; // (The state of IkB)
	in >> _GR; // (The state of reporter genes)
	in >> _c1rrChemTNF; // NF-kB independent rate of TNF/chemokine mRNA synthesis  
	in >> _c1rChem;
	in >> _c1rTNF;
	in >> _chemt; // generic chemokine transcript
	in >> _chem; // intracellular generic chemokine protein
	in >> _TNFt; // TNF transcript
	in >> _TNF; // intracellular TNF 
	in >> _ACTt; // transcript of macrophage activating molecules
	in >> _ACT; // mac-activation molecules
	in >> _normalizedACT;
	in >> _IAPt; // transcript of IAP (inhibitor of apoptosis)
	in >> _IAP; 
	in >> _normalizedIAP;
    
    size_t vectorSize;
    in >> vectorSize;

    _initvector.resize(vectorSize, 0.0);
    _k1vector.resize(vectorSize, 0.0);
    _k2vector.resize(vectorSize, 0.0);
    _k3vector.resize(vectorSize, 0.0);
    _k4vector.resize(vectorSize, 0.0);
    _switchvector.resize(vectorSize, 0.0);

    for (size_t kk = 0; kk < vectorSize; kk++) 
    {
        in >> _initvector[kk];
        in >> _k1vector[kk];
        in >> _k2vector[kk];
        in >> _k3vector[kk];
        in >> _k4vector[kk];
        in >> _switchvector[kk];
    }
    
	if (!Serialization::readFooter(in, Agent::_ClassName))
	{
		exit(1);
	}
}
