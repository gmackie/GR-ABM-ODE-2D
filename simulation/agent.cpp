/*
 * agent.cpp
 *
 *  Created on: 02-nov-2009
 *      Author: M. El-Kebir
 */

#include "agent.h"
#include "grgrid.h"

using namespace std;

auto_ptr<ODESolvers::Stepper> Agent::stepper;
auto_ptr<LungFunc> Agent::deriv;

unsigned long Agent::_nextID = 0;

// Needed for deserializing the model state.
// Avoids the calls to the random number generator in the normal constructor, allowing the random number generator
// to remain in synch after deserialization.
Agent::Agent()
  :
#define P(type, name, ival, desc) \
      _##name (),
  AGENT_PROPS
#undef P
  _lasttimestep(-1.0),
  _initvector()
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
            )
  :
#define P(type, name, ival, desc) \
      _##name (ival) ,
  AGENT_PROPS
#undef P
  _lasttimestep(_PARAM(_timestepMolecular))
{
  _id = createID();
  _birthTime = (birthtime);
  _deathTime = (deathtime);
  _pos = Pos(row, col);

  // TNFR components
  _surfTNFR1 = (g_Rand.getReal(meanTNFR1 - stdTNFR1, meanTNFR1 + stdTNFR1));
  _surfTNFR2 = (g_Rand.getReal(meanTNFR2 - stdTNFR2, meanTNFR2 + stdTNFR2));
  //_surfTNFR1 = (g_Rand.getLogNormal(meanTNFR1),_PARAM(stdTNFR1)));
  //_surfTNFR2 = (g_Rand.getLogNormal(meanTNFR2),_PARAM(stdTNFR2)));
  _vTNFR1 = (_surfTNFR1 * _PARAM(_kT1));
  _vTNFR2 = (_surfTNFR2 * _PARAM(_kT2));
  _kSynth = (kSynth);
  _kTACE = (kTACE);

  // IL10 components
  _surfIL10R = (g_Rand.getReal(iIL10R - stdIL10R, iIL10R + stdIL10R));
  _vIL10R = (_surfIL10R * _PARAM(_IkT));
  _meanIL10R = (iIL10R);

  // NF-kB signaling pathway components
  _NFkB_IkB = ((_PARAM(_meanNFkB) > 0 ) ? g_Rand.getLogNormal(_PARAM(_meanNFkB),0.64872*_PARAM(_meanNFkB)) : 0);
}

Agent::~Agent()
{
}

void Agent::solveMolecularScale(GrGrid& grid, double t, double dt, ODESolvers::ODEMethod method)
{
  LungFunc& fn = *getDerivFunc();
  if(_initvector.size() != fn.dim()) _initvector.resize(fn.dim());
  writeValarrayFromMembers(grid, _initvector);  //With some pointer magic, should be able to remove this...
  static valarray<double> error(fn.dim()); //For adaptive only, not needed atm
  static LungFunc::Params_t params;
  params.agent = this;
  params.grid = &grid;
  size_t nsubsteps = (_PARAM(_NFkBdynamics) && NFkBCapable()) ? _PARAM(_NFkBTimeCoefficient) : 1;
  ODESolvers::Stepper* stepper = getStepper(method);
  for(size_t i=0; i<nsubsteps; i++)
    stepper->step(_initvector, fn, t+i*dt/nsubsteps, dt/nsubsteps, _lasttimestep, error, (void*)&params);
  writeMembersFromValarray(grid, _initvector);
}

void Agent::writeValarrayFromMembers(GrGrid& grid, valarray<double>& inputVector)
{
  // Defining the initial vector from the agent private variables
  // and getting soluble values from the grid
  int vectorsize = inputVector.size();

  if (vectorsize == 3)
    {
      // sIL10
      inputVector[0] = grid.il10(_pos);
//        inputVector[0] = grid.il10(_pos) / (NAV * VOL);
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
      inputVector[8] = grid.TNF(_pos);
//        inputVector[8] = grid.TNF(_pos) / (NAV * VOL);
      // shedTNFR2
      inputVector[9] = grid.shedTNFR2(_pos);
//        inputVector[9] = grid.shedTNFR2(_pos) / (NAV * VOL);

      if (vectorsize == 13)
        {
          // sIL10
          inputVector[10] = grid.il10(_pos);
//            inputVector[10] = grid.il10(_pos) / (NAV * VOL);
          // surfIL10R
          inputVector[11] = _surfIL10R;
          // surfBoundIL10R
          inputVector[12] = _surfBoundIL10R;
        }

      if (vectorsize == 36 || vectorsize == 39)
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
          inputVector[34] = _IAPt;
          inputVector[35] = _IAP;
#if 0
          inputVector[36] = _normalizedACT;
          inputVector[37] = _normalizedIAP;
#endif

          if (vectorsize == 39)
            {
              // sIL10
              inputVector[36] = grid.il10(_pos);
//                inputVector[36] = grid.il10(_pos) / (NAV * VOL);
              // surfIL10R
              inputVector[37] = _surfIL10R;
              // surfBoundIL10R
              inputVector[38] = _surfBoundIL10R;
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
      grid.setil10(_pos, (inputVector[0]));
//        grid.setil10(_pos, (NAV * VOL * inputVector[0]));
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
      grid.setTNF(_pos, (inputVector[8]));
//        grid.setTNF(_pos, (NAV * VOL * inputVector[8]));
      // shedTNFR2
      grid.setshedTNFR2(_pos, (inputVector[9]));
//        grid.setshedTNFR2(_pos, (NAV * VOL * inputVector[9]));

      if (vectorsize == 13)
        {
          // sIL10
          grid.setil10(_pos, (inputVector[10]));
//            grid.setil10(_pos, (NAV * VOL * inputVector[10]));
          // surfIL10R
          _surfIL10R = inputVector[11];
          // surfBoundIL10R
          _surfBoundIL10R = inputVector[12];
        }

      if (vectorsize == 36 || vectorsize == 39)
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
          _IAPt = inputVector[34];
          _IAP = inputVector[35];
#if 0
          _normalizedACT = inputVector[36];
          _normalizedIAP = inputVector[37];
#else
          _normalizedACT = _ACT*_PARAM(_c3rACT);
          _normalizedIAP = _IAP*_PARAM(_c3rIAP);
#endif

          grid.incCCL2(_pos, (_PARAM(_e3Chem) * inputVector[29] * _PARAM(_timestepMolecular)));
          grid.incCCL5(_pos, (_PARAM(_e3Chem) * inputVector[29] * _PARAM(_timestepMolecular)));
          grid.incCXCL9(_pos, (2 * _PARAM(_e3Chem) * inputVector[29] * _PARAM(_timestepMolecular)));

          if (vectorsize == 39)
            {
              // sIL10
              grid.setil10(_pos, (inputVector[36]));
//                grid.setil10(_pos, (NAV * VOL * inputVector[36]));
              // surfIL10R
              _surfIL10R = inputVector[37];
              // surfBoundIL10R
              _surfBoundIL10R = inputVector[38];
            }
        }
    }
}

void Agent::checkTolerance(valarray<double>& veccheck)
{
  double intpart, TempStoreLarge, TempStorePower;
  int intpartStore;

  // Checks sig figs and round correctly to that number
  for (int i = 0; i < (int)veccheck.size(); i++)
    {
      if (veccheck[i] > 0)
        {
          TempStorePower = floor(log10(veccheck[i]));
          TempStoreLarge = (veccheck[i] * (ABS_TOL/(pow(10,TempStorePower))));
          modf(TempStoreLarge, &intpart);
          intpartStore = (int)intpart;
//            if (fracpart >= 0.5)
//            {
//                intpartStore += (int)(ceil(fracpart));
//            }
          TempStoreLarge = (intpartStore / (ABS_TOL/(pow(10,TempStorePower))));
          veccheck[i] = TempStoreLarge;
        }
    }
}

bool Agent::intCompareGT(const double param1, const double param2)
{
  // Compares two values based on the number of sig figs we hold in gr.h (ABS_TOL)
  // If the values are not within 2 orders of magnitude we do not convert to ints
  // since it should not matter
  // COMPARES PARAM1 > PARAM2 and returns bool based on this evaluation

  double intpart1, intpart2, Store1, Store2, StorePower;
  int intpart1Store, intpart2Store;

  bool result = 0;

  if (fabs(floor(log10(param1)) - floor(log10(param2))) < 2  )
    {
      if (floor(log10(param1)) < floor(log10(param2)))
        {
          StorePower = floor(log10(param1));
        }
      else
        {
          StorePower = floor(log10(param2));
        }

      Store1 = (param1 * (ABS_TOL/(pow(10,StorePower))));
      Store2 = (param2 * (ABS_TOL/(pow(10,StorePower))));

      modf(Store1, &intpart1);
      modf(Store2, &intpart2);

      intpart1Store = (int)intpart1;
      intpart2Store = (int)intpart2;

      if (intpart1Store > intpart2Store)
        {
          result = 1;
        }

//        std::cout << param1 << "   " << param2 << std::endl;
//        std::cout << intpart1Store << "   " << intpart2Store << std::endl;

    }
  else
    {
      if (param1 > param2)
        {
          result = 1;
        }

//        std::cout << param1 << "   " << param2 << std::endl;
    }

//    std::cout << result << std::endl;

  return result;
}

#if 0

void Agent::derivativeTNF(const valarray<double>& vecread, valarray<double>& vecwrite, double dt, GrGrid& grid)
{
  assert(vecread.size() == 10); // Make sure the valarray length is set correctly

  double koff1 = _PARAM(_kOn1) * _PARAM(_KD1);
  double koff2 = _PARAM(_kOn2) * _PARAM(_KD2);
  double il10 = grid.il10(_pos) /(NAV * VOL);
  double IkmRNA;

  // modulate TNF parameters based on equilibrium IL10 receptor equations since the il10 molecular equations are not active
  double eqsurfBoundIL10R = (il10 * _meanIL10R) / (_PARAM(_IkD) + il10);

  if (_kmRNA > 0)
    {
      IkmRNA = _kmRNA * ((_kSynth/_kmRNA)+ ((1.0 - (_kSynth/_kmRNA))/(1.0 + exp((eqsurfBoundIL10R - _PARAM(_LinkRNAGamma))/_PARAM(_LinkRNADelta)))));
    }
  else
    {
      IkmRNA = 0.0;
    }

  // TNF Ordinary Differential Equations
  // mTNFRNA
  vecwrite[0] = (IkmRNA - _PARAM(_kTrans) * vecread[0]) * dt;
  // mTNF
  vecwrite[1] = (_PARAM(_kTrans) * vecread[0] - _kTACE * vecread[1]) * dt;
  // surfTNFR1
  vecwrite[2] = (_vTNFR1 - _PARAM(_kOn1) * vecread[8] * vecread[2] + koff1 * vecread[4] - _PARAM(_kT1) * vecread[2] + _PARAM(_kRec1) * vecread[6]) * dt;
  // surfTNFR2
  vecwrite[3] = (_vTNFR2 - _PARAM(_kOn2) * vecread[8] * vecread[3] + koff2 * vecread[5] - _PARAM(_kT2) * vecread[3] + _PARAM(_kRec2) * vecread[7]) * dt;
  // surfBoundTNFR1
  vecwrite[4] = (_PARAM(_kOn1) * vecread[8] * vecread[2] - koff1 * vecread[4] - _PARAM(_kInt1) * vecread[4]) * dt;
  // surfBoundTNFR2
  vecwrite[5] = (_PARAM(_kOn2) * vecread[8] * vecread[3] - koff2 * vecread[5] - _PARAM(_kInt2) * vecread[5] - _PARAM(_kShed) * vecread[5]) * dt;
  // intBoundTNFR1
  vecwrite[6] = (_PARAM(_kInt1) * vecread[4] - _PARAM(_kDeg1) * vecread[6] - _PARAM(_kRec1) * vecread[6]) * dt;
  // intBoundTNFR2
  vecwrite[7] = (_PARAM(_kInt2) * vecread[5] - _PARAM(_kDeg2) * vecread[7] - _PARAM(_kRec2) * vecread[7]) * dt;
  // sTNF
  vecwrite[8] = ((DENSITY/NAV) * (_kTACE * vecread[1] - _PARAM(_kOn1) * vecread[8] * vecread[2] + koff1 * vecread[4] - _PARAM(_kOn2) * vecread[8] * vecread[3] + koff2 * vecread[5])) * dt;
  // shedTNFR2
  vecwrite[9] = ((DENSITY/NAV) * _PARAM(_kShed) * vecread[5]) * dt;
}

void Agent::solveTNF(GrGrid& grid, double dt)
{
  double koff1 = _PARAM(_kOn1) * _PARAM(_KD1);
  double koff2 = _PARAM(_kOn2) * _PARAM(_KD2);

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
  eqsurfBoundIL10R = (il10 * _meanIL10R) / (_PARAM(_IkD) + il10);

  if (_kmRNA > 0)
    {
      IkmRNA = _kmRNA * ((_kSynth/_kmRNA)+ ((1.0 - (_kSynth/_kmRNA))/(1.0 + exp((eqsurfBoundIL10R - _PARAM(_LinkRNAGamma))/_PARAM(_LinkRNADelta)))));
    }
  else
    {
      IkmRNA = 0.0;
    }

  // end of equilibrium calculations

  dmTNFRNA = (IkmRNA - _PARAM(_kTrans) * _mTNFRNA) * dt;
  dmTNF = (_PARAM(_kTrans) * _mTNFRNA - _kTACE * _mTNF) * dt;
  dsurfTNFR1 = (_vTNFR1 - _PARAM(_kOn1) * tnf * _surfTNFR1 + koff1 * _surfBoundTNFR1 - _PARAM(_kT1) * _surfTNFR1 + _PARAM(_kRec1) * _intBoundTNFR1) * dt;
  dsurfTNFR2 = (_vTNFR2 - _PARAM(_kOn2) * tnf * _surfTNFR2 + koff2 * _surfBoundTNFR2 - _PARAM(_kT2) * _surfTNFR2 + _PARAM(_kRec2) * _intBoundTNFR2) * dt;
  dsurfBoundTNFR1 = (_PARAM(_kOn1) * tnf * _surfTNFR1 - koff1 * _surfBoundTNFR1 - _PARAM(_kInt1) * _surfBoundTNFR1) * dt;
  dsurfBoundTNFR2 = (_PARAM(_kOn2) * tnf * _surfTNFR2 - koff2 * _surfBoundTNFR2 - _PARAM(_kInt2) * _surfBoundTNFR2 - _PARAM(_kShed) * _surfBoundTNFR2) * dt;
  dintBoundTNFR1 = (_PARAM(_kInt1) * _surfBoundTNFR1 - _PARAM(_kDeg1) * _intBoundTNFR1 - _PARAM(_kRec1) * _intBoundTNFR1) * dt;
  dintBoundTNFR2 = (_PARAM(_kInt2) * _surfBoundTNFR2 - _PARAM(_kDeg2) * _intBoundTNFR2 - _PARAM(_kRec2) * _intBoundTNFR2) * dt;
  dsTNF = ((DENSITY/NAV) * (_kTACE * _mTNF - _PARAM(_kOn1) * tnf * _surfTNFR1 + koff1 * _surfBoundTNFR1 - _PARAM(_kOn2) * tnf * _surfTNFR2 + koff2 * _surfBoundTNFR2)) * dt;
  dshedTNFR2 = ((DENSITY/NAV) * _PARAM(_kShed) * _surfBoundTNFR2) * dt;

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

  assert(vecread.size() == 3); // Make sure the valarray length is set correctly

  double Ikoff = _PARAM(_IkOn) * _PARAM(_IkD);


  // IL10 Ordinary Differential Equations
  // sIL10
  vecwrite[0] = (((DENSITY/NAV) * _kISynth) + ((DENSITY/NAV) * (Ikoff * vecread[2] - _PARAM(_IkOn) * vecread[1] * vecread[0]))) * dt;
  // surfIL10R
  vecwrite[1] = (_vIL10R - _PARAM(_IkOn) * vecread[1] * vecread[0] + _PARAM(_IkOff) * vecread[2] - _PARAM(_IkT) * vecread[1]) * dt;
  // surfBoundIL10R
  vecwrite[2] = (_PARAM(_IkOn) * vecread[1] * vecread[0] - _PARAM(_IkOff) * vecread[2] - _PARAM(_IkInt) * vecread[2]) * dt;

}

void Agent::solveIL10(GrGrid& grid, double dt)
{
  double il10 = grid.il10(_pos) / (NAV * VOL);
  double Ikoff = _PARAM(_IkOn) * _PARAM(_IkD);

  double dsIL10;
  double dsurfIL10R;
  double dsurfBoundIL10R;

  // IL10 differential equations
  dsIL10 = ((DENSITY/NAV) * _kISynth + ((DENSITY/NAV) * (Ikoff * _surfBoundIL10R - _PARAM(_IkOn) * _surfIL10R * il10))) * dt;
  dsurfIL10R = (_vIL10R - _PARAM(_IkOn) * _surfIL10R * il10 + Ikoff * _surfBoundIL10R - _PARAM(_IkT) * _surfIL10R) * dt;
  dsurfBoundIL10R = (_PARAM(_IkOn) * _surfIL10R * il10 - Ikoff * _surfBoundIL10R - _PARAM(_IkInt) * _surfBoundIL10R) * dt;
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

void Agent::derivativeTNFandIL10(const valarray<double>& vecread, valarray<double>& vecwrite, double dt, GrGrid& grid)
{
  assert(vecread.size() == 13); // Make sure the valarray length is set correctly

  double koff1 = _PARAM(_kOn1) * _PARAM(_KD1);
  double koff2 = _PARAM(_kOn2) * _PARAM(_KD2);
  double Ikoff = _PARAM(_IkOn) * _PARAM(_IkD);
  double IkmRNA;

  // solving for TNF parameters that depend on IL10

  if (_kmRNA > 0)
    {
      IkmRNA = _kmRNA * ((_kSynth/_kmRNA) + ((1.0 - (_kSynth/_kmRNA))/(1.0 + exp((vecread[12] - _PARAM(_LinkRNAGamma))/_PARAM(_LinkRNADelta)))));
    }
  else
    {
      IkmRNA = 0.0;
    }

  // TNF Ordinary Differential Equations
  // mTNFRNA
  vecwrite[0] = (IkmRNA - _PARAM(_kTrans) * vecread[0]) * dt;
  // mTNF
  vecwrite[1] = (_PARAM(_kTrans) * vecread[0] - _kTACE * vecread[1]) * dt;
  // surfTNFR1
  vecwrite[2] = (_vTNFR1 - _PARAM(_kOn1) * vecread[8] * vecread[2] + koff1 * vecread[4] - _PARAM(_kT1) * vecread[2] + _PARAM(_kRec1) * vecread[6]) * dt;
  // surfTNFR2
  vecwrite[3] = (_vTNFR2 - _PARAM(_kOn2) * vecread[8] * vecread[3] + koff2 * vecread[5] - _PARAM(_kT2) * vecread[3] + _PARAM(_kRec2) * vecread[7]) * dt;
  // surfBoundTNFR1
  vecwrite[4] = (_PARAM(_kOn1) * vecread[8] * vecread[2] - koff1 * vecread[4] - _PARAM(_kInt1) * vecread[4]) * dt;
  // surfBoundTNFR2
  vecwrite[5] = (_PARAM(_kOn2) * vecread[8] * vecread[3] - koff2 * vecread[5] - _PARAM(_kInt2) * vecread[5] - _PARAM(_kShed) * vecread[5]) * dt;
  // intBoundTNFR1
  vecwrite[6] = (_PARAM(_kInt1) * vecread[4] - _PARAM(_kDeg1) * vecread[6] - _PARAM(_kRec1) * vecread[6]) * dt;
  // intBoundTNFR2
  vecwrite[7] = (_PARAM(_kInt2) * vecread[5] - _PARAM(_kDeg2) * vecread[7] - _PARAM(_kRec2) * vecread[7]) * dt;
  // sTNF
  vecwrite[8] = ((DENSITY/NAV) * (_kTACE * vecread[1] - _PARAM(_kOn1) * vecread[8] * vecread[2] + koff1 * vecread[4] - _PARAM(_kOn2) * vecread[8] * vecread[3] + koff2 * vecread[5])) * dt;
  // shedTNFR2
  vecwrite[9] = ((DENSITY/NAV) * _PARAM(_kShed) * vecread[5]) * dt;

  // IL10 Ordinary Differential Equations
  // sIL10
  vecwrite[10] = (((DENSITY/NAV) * _kISynth) + ((DENSITY/NAV) * (Ikoff * vecread[12] - _PARAM(_IkOn) * vecread[11] * vecread[10]))) * dt;
  // surfIL10R
  vecwrite[11] = (_vIL10R - _PARAM(_IkOn) * vecread[11] * vecread[10] + Ikoff * vecread[12] - _PARAM(_IkT) * vecread[11]) * dt;
  // surfBoundIL10R
  vecwrite[12] = (_PARAM(_IkOn) * vecread[11] * vecread[10] - Ikoff * vecread[12] - _PARAM(_IkInt) * vecread[12]) * dt;

}

void Agent::solveTNFandIL10(GrGrid& grid, double dt)
{
  double koff1 = _PARAM(_kOn1) * _PARAM(_KD1);
  double koff2 = _PARAM(_kOn2) * _PARAM(_KD2);
  double Ikoff = _PARAM(_IkOn) * _PARAM(_IkD);

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
      IkmRNA = _kmRNA * ((_kSynth/_kmRNA) + ((1.0 - (_kSynth/_kmRNA))/(1.0 + exp((_surfBoundIL10R - _PARAM(_LinkRNAGamma))/_PARAM(_LinkRNADelta)))));
    }
  else
    {
      IkmRNA = 0.0;
    }

  // end of TNF and IL10 linking

  // TNF differential equations
  dmTNFRNA = (IkmRNA - _PARAM(_kTrans) * _mTNFRNA) * dt;
  dmTNF = (_PARAM(_kTrans) * _mTNFRNA - _kTACE * _mTNF) * dt;
  dsurfTNFR1 = (_vTNFR1 - _PARAM(_kOn1) * tnf * _surfTNFR1 + koff1 * _surfBoundTNFR1 - _PARAM(_kT1) * _surfTNFR1 + _PARAM(_kRec1) * _intBoundTNFR1) * dt;
  dsurfTNFR2 = (_vTNFR2 - _PARAM(_kOn2) * tnf * _surfTNFR2 + koff2 * _surfBoundTNFR2 - _PARAM(_kT2) * _surfTNFR2 + _PARAM(_kRec2) * _intBoundTNFR2) * dt;
  dsurfBoundTNFR1 = (_PARAM(_kOn1) * tnf * _surfTNFR1 - koff1 * _surfBoundTNFR1 - _PARAM(_kInt1) * _surfBoundTNFR1) * dt;
  dsurfBoundTNFR2 = (_PARAM(_kOn2) * tnf * _surfTNFR2 - koff2 * _surfBoundTNFR2 - _PARAM(_kInt2) * _surfBoundTNFR2 - _PARAM(_kShed) * _surfBoundTNFR2) * dt;
  dintBoundTNFR1 = (_PARAM(_kInt1) * _surfBoundTNFR1 - _PARAM(_kDeg1) * _intBoundTNFR1 - _PARAM(_kRec1) * _intBoundTNFR1) * dt;
  dintBoundTNFR2 = (_PARAM(_kInt2) * _surfBoundTNFR2 - _PARAM(_kDeg2) * _intBoundTNFR2 - _PARAM(_kRec2) * _intBoundTNFR2) * dt;
  dsTNF = ((DENSITY/NAV) * (_kTACE * _mTNF - _PARAM(_kOn1) * tnf * _surfTNFR1 + koff1 * _surfBoundTNFR1 - _PARAM(_kOn2) * tnf * _surfTNFR2 + koff2 * _surfBoundTNFR2)) * dt;
  dshedTNFR2 = ((DENSITY/NAV) * _PARAM(_kShed) * _surfBoundTNFR2) * dt;
  // end of TNF differential equations

  // IL10 differential equations
  dsIL10 = ((DENSITY/NAV) * _kISynth + ((DENSITY/NAV) * (Ikoff * _surfBoundIL10R - _PARAM(_IkOn) * _surfIL10R * il10))) * dt;
  dsurfIL10R = (_vIL10R - _PARAM(_IkOn) * _surfIL10R * il10 + Ikoff * _surfBoundIL10R - _PARAM(_IkT) * _surfIL10R) * dt;
  dsurfBoundIL10R = (_PARAM(_IkOn) * _surfIL10R * il10 - Ikoff * _surfBoundIL10R - _PARAM(_IkInt) * _surfBoundIL10R) * dt;
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

void Agent::derivativeTNFandNFKB(const valarray<double>& vecread, valarray<double>& vecwrite, double dt, GrGrid& grid)
{
  assert(vecread.size() == 38);

  double koff1 = _PARAM(_kOn1) * _PARAM(_KD1);
  double koff2 = _PARAM(_kOn2) * _PARAM(_KD2);
  double il10 = grid.il10(_pos) /(NAV * VOL);
  double eqsurfBoundIL10R;
  double IkmRNA;

  // modulate TNF parameters based on equilibrium IL10 receptor equations since the il10 molecular equations are not active
  eqsurfBoundIL10R = (il10 * _surfIL10R) / (_PARAM(_IkD) + il10);

  if (_kmRNA > 0)
    {
      IkmRNA = _kmRNA * ((_kSynth/_kmRNA)+ ((1.0 - (_kSynth/_kmRNA))/(1.0 + exp((eqsurfBoundIL10R - _PARAM(_LinkRNAGamma))/_PARAM(_LinkRNADelta)))));
    }
  else
    {
      IkmRNA = 0.0;
    }

  // end of equilibrium calculations


  // CURRENTLY NOT LINKED WITH IL10 SINCE NFKB DYNAMICS HAVE NOT BEEN LOOKED AT

  // TNF Ordinary Differential Equations

  vecwrite[0] = 0 * dt; // CURRENTLY NOT LINKED WITH IL10 SINCE NFKB DYNAMICS HAVE NOT BEEN LOOKED AT
  vecwrite[1] = (_PARAM(_e3TNF)*_TNF - _kTACE * vecread[1]) * dt;
  vecwrite[2] = (_vTNFR1 - _PARAM(_kOn1) * vecread[8] * vecread[2] + (koff1+KDEG) * vecread[4] - _PARAM(_kT1) * vecread[2] + _PARAM(_kRec1) * vecread[6]) * dt;
  vecwrite[3] = (_vTNFR2 - _PARAM(_kOn2) * vecread[8] * vecread[3] + (koff2+KDEG) * vecread[5] - _PARAM(_kT2) * vecread[3] + _PARAM(_kRec2) * vecread[7]) * dt;
  vecwrite[4] = (_PARAM(_kOn1) * vecread[8] * vecread[2] - (koff1+KDEG) * vecread[4] - _PARAM(_kInt1) * vecread[4]) * dt;
  vecwrite[5] = (_PARAM(_kOn2) * vecread[8] * vecread[3] - (koff2+KDEG) * vecread[5] - _PARAM(_kInt2) * vecread[5] - _PARAM(_kShed) * vecread[5]) * dt;
  vecwrite[6] = (_PARAM(_kInt1) * vecread[4] - _PARAM(_kDeg1) * vecread[6] - _PARAM(_kRec1) * vecread[6]) * dt;
  vecwrite[7] = (_PARAM(_kInt2) * vecread[5] - _PARAM(_kDeg2) * vecread[7] - _PARAM(_kRec2) * vecread[7]) * dt;
  vecwrite[8] = ((DENSITY/NAV) * (_kTACE * vecread[1] - _PARAM(_kOn1) * vecread[8] * vecread[2] + koff1 * vecread[4] - _PARAM(_kOn2) * vecread[8] * vecread[3] + koff2 * vecread[5])) * dt;
  vecwrite[9] = ((DENSITY/NAV) * _PARAM(_kShed) * vecread[5]) * dt;
  // NF-kB dynamics model equations
  vecwrite[10] = (_PARAM(_ka)*vecread[4]*(_PARAM(_KN)-vecread[10])*_PARAM(_kA20)/(_PARAM(_kA20)+vecread[18])-_PARAM(_ki)*vecread[10]) * dt;
  vecwrite[11] = (_PARAM(_k4)*(_PARAM(_KNN)-vecread[11]-vecread[12]-vecread[13])-_PARAM(_k1)*(vecread[10]*vecread[10])*vecread[11]) * dt;
  vecwrite[12] = (_PARAM(_k1)*(vecread[10]*vecread[10])*vecread[11]-_PARAM(_k3)*vecread[12]*(_PARAM(_k2)+vecread[18])/_PARAM(_k2)) * dt;
  vecwrite[13] = (_PARAM(_k3)*vecread[12]*(_PARAM(_k2)+vecread[18])/_PARAM(_k2)-_PARAM(_k4)*vecread[13]) * dt;
  vecwrite[14] = (_PARAM(_a2)*vecread[12]*vecread[20]-_PARAM(_tp)*vecread[14]) * dt;
  vecwrite[15] = (_PARAM(_a3)*vecread[12]*vecread[23]-_PARAM(_tp)*vecread[15]) * dt;
  vecwrite[16] = (_PARAM(_c6a)*vecread[23]-_PARAM(_a1)*vecread[16]*vecread[20]+_PARAM(_tp)*vecread[15]-_PARAM(_i1)*vecread[16]) * dt;
  vecwrite[17] = (_PARAM(_i1)*vecread[16]-_PARAM(_a1)*KV*vecread[21]*vecread[17]) * dt;
  vecwrite[18] = (_PARAM(_c4)*vecread[19]-_PARAM(_c5)*vecread[18]) * dt;
  vecwrite[19] = (_PARAM(_c1)*vecread[25]-_PARAM(_c3)*vecread[19]) * dt;
  vecwrite[20] = (-_PARAM(_a2)*vecread[12]*vecread[20]-_PARAM(_a1)*vecread[20]*vecread[16]+_PARAM(_c4)*vecread[22]-_PARAM(_c5a)*vecread[20]-_PARAM(_i1a)*vecread[20]+_PARAM(_e1a)*vecread[21]) * dt;
  vecwrite[21] = (-_PARAM(_a1)*KV*vecread[21]*vecread[17]+_PARAM(_i1a)*vecread[20]-_PARAM(_e1a)*vecread[21]) * dt;
  vecwrite[22] = (_PARAM(_c1)*vecread[26]-_PARAM(_c3)*vecread[22]) * dt;
  vecwrite[23] = (_PARAM(_a1)*vecread[20]*vecread[16]-_PARAM(_c6a)*vecread[23]-_PARAM(_a3)*vecread[12]*vecread[23]+_PARAM(_e2a)*vecread[24]) * dt;
  vecwrite[24] = (_PARAM(_a1)*KV*vecread[21]*vecread[17]-_PARAM(_e2a)*vecread[24]) * dt;
  vecwrite[25] = (_PARAM(_q1)*vecread[17]*(2-vecread[25])-_PARAM(_q2)*vecread[21]*vecread[25]) * dt;
  vecwrite[26] = (_PARAM(_q1)*vecread[17]*(2-vecread[26])-_PARAM(_q2)*vecread[21]*vecread[26]) * dt;
  vecwrite[27] = (_PARAM(_q1r)*vecread[17]*(2-vecread[27])-(_PARAM(_q2rr)+_PARAM(_q2r)*vecread[21])*vecread[27]) * dt;
  vecwrite[28] = (_c1rrChemTNF+_c1rChem*vecread[27]-_PARAM(_c3rChem)*vecread[28]) * dt;
  vecwrite[29] = (_PARAM(_c4Chem)*vecread[28]-_PARAM(_c5Chem)*vecread[29]-_PARAM(_e3Chem)*vecread[29]) * dt;
  vecwrite[30] = (_c1rrChemTNF+_c1rTNF*vecread[27]-_PARAM(_c3rTNF)*vecread[30]) * dt;
  vecwrite[31] = (_PARAM(_c4TNF)*vecread[30]-_PARAM(_c5TNF)*vecread[31]-_PARAM(_e3TNF)*vecread[31]) * dt;
  vecwrite[32] = (_PARAM(_c1rrACT)+_PARAM(_c1r)*vecread[27]-_PARAM(_c3rACT)*vecread[32]) * dt;
  vecwrite[33] = (_PARAM(_c4ACT)*vecread[32]-_PARAM(_c5ACT)*vecread[33]) * dt;
  vecwrite[34] = (_PARAM(_c1rrIAP)+_PARAM(_c1r)*vecread[27]-_PARAM(_c3rIAP)*vecread[34]) * dt;
  vecwrite[35] = (_PARAM(_c4IAP)*vecread[34]-_PARAM(_c5IAP)*vecread[35]) * dt;
  vecwrite[36] = vecread[33]*_PARAM(_c3rACT);
  vecwrite[37] = vecread[35]*_PARAM(_c3rIAP);

}

void Agent::solveNFkBandTNF(GrGrid& grid, double dt)
{
  double koff1 = _PARAM(_kOn1) * _PARAM(_KD1);
  double koff2 = _PARAM(_kOn2) * _PARAM(_KD2);

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
  eqsurfBoundIL10R = (il10 * _surfIL10R) / (_PARAM(_IkD) + il10);

  if (_kmRNA > 0)
    {
      IkmRNA = _kmRNA * ((_kSynth/_kmRNA)+ ((1.0 - (_kSynth/_kmRNA))/(1.0 + exp((eqsurfBoundIL10R - _PARAM(_LinkRNAGamma))/_PARAM(_LinkRNADelta)))));
    }
  else
    {
      IkmRNA = 0.0;
    }

  // end of equilibrium calculations


  // CURRENTLY NOT LINKED WITH IL10 SINCE NFKB DYNAMICS HAVE NOT BEEN LOOKED AT

  // tnf dynamics model equations
  dmTNF = (_PARAM(_e3TNF)*_TNF - _kTACE * _mTNF) * dt;
  dsurfTNFR1 = (_vTNFR1 - _PARAM(_kOn1) * tnf * _surfTNFR1 + (koff1+KDEG) * _surfBoundTNFR1 - _PARAM(_kT1) * _surfTNFR1 + _PARAM(_kRec1) * _intBoundTNFR1) * dt;
  dsurfTNFR2 = (_vTNFR2 - _PARAM(_kOn2) * tnf * _surfTNFR2 + (koff2+KDEG) * _surfBoundTNFR2 - _PARAM(_kT2) * _surfTNFR2 + _PARAM(_kRec2) * _intBoundTNFR2) * dt;
  dsurfBoundTNFR1 = (_PARAM(_kOn1) * tnf * _surfTNFR1 - (koff1+KDEG) * _surfBoundTNFR1 - _PARAM(_kInt1) * _surfBoundTNFR1) * dt;
  dsurfBoundTNFR2 = (_PARAM(_kOn2) * tnf * _surfTNFR2 - (koff2+KDEG) * _surfBoundTNFR2 - _PARAM(_kInt2) * _surfBoundTNFR2 - _PARAM(_kShed) * _surfBoundTNFR2) * dt;
  dintBoundTNFR1 = (_PARAM(_kInt1) * _surfBoundTNFR1 - _PARAM(_kDeg1) * _intBoundTNFR1 - _PARAM(_kRec1) * _intBoundTNFR1) * dt;
  dintBoundTNFR2 = (_PARAM(_kInt2) * _surfBoundTNFR2 - _PARAM(_kDeg2) * _intBoundTNFR2 - _PARAM(_kRec2) * _intBoundTNFR2) * dt;
  dsTNF = ((DENSITY/NAV) * (_kTACE * _mTNF - _PARAM(_kOn1) * tnf * _surfTNFR1 + koff1 * _surfBoundTNFR1 - _PARAM(_kOn2) * tnf * _surfTNFR2 + koff2 * _surfBoundTNFR2)) * dt;
  dshedTNFR2 = ((DENSITY/NAV) * _PARAM(_kShed) * _surfBoundTNFR2) * dt;

  // NF-kB dynamics model equations
  dIKKKa = (_PARAM(_ka)*_surfBoundTNFR1*(_PARAM(_KN)-_IKKKa)*_PARAM(_kA20)/(_PARAM(_kA20)+_A20)-_PARAM(_ki)*_IKKKa) * dt;
  dIKKn = (_PARAM(_k4)*(_PARAM(_KNN)-_IKKn-_IKKa-_IKKi)-_PARAM(_k1)*(_IKKKa*_IKKKa)*_IKKn) * dt;
  dIKKa = (_PARAM(_k1)*(_IKKKa*_IKKKa)*_IKKn-_PARAM(_k3)*_IKKa*(_PARAM(_k2)+_A20)/_PARAM(_k2)) * dt;
  dIKKi = (_PARAM(_k3)*_IKKa*(_PARAM(_k2)+_A20)/_PARAM(_k2)-_PARAM(_k4)*_IKKi) * dt;
  dIkBp = (_PARAM(_a2)*_IKKa*_IkB-_PARAM(_tp)*_IkBp) * dt;
  dNFkB_IkBp = (_PARAM(_a3)*_IKKa*_NFkB_IkB-_PARAM(_tp)*_NFkB_IkBp) * dt;
  dNFkBc = (_PARAM(_c6a)*_NFkB_IkB-_PARAM(_a1)*_NFkBc*_IkB+_PARAM(_tp)*_NFkB_IkBp-_PARAM(_i1)*_NFkBc) * dt;
  dNFkBn = (_PARAM(_i1)*_NFkBc-_PARAM(_a1)*KV*_IkBn*_NFkBn) * dt;
  dA20 = (_PARAM(_c4)*_A20t-_PARAM(_c5)*_A20) * dt;
  dA20t = (_PARAM(_c1)*_GA20-_PARAM(_c3)*_A20t) * dt;
  dIkB = (-_PARAM(_a2)*_IKKa*_IkB-_PARAM(_a1)*_IkB*_NFkBc+_PARAM(_c4)*_IkBt-_PARAM(_c5a)*_IkB-_PARAM(_i1a)*_IkB+_PARAM(_e1a)*_IkBn) * dt;
  dIkBn = (-_PARAM(_a1)*KV*_IkBn*_NFkBn+_PARAM(_i1a)*_IkB-_PARAM(_e1a)*_IkBn) * dt;
  dIkBt = (_PARAM(_c1)*_GIkB-_PARAM(_c3)*_IkBt) * dt;
  dNFkB_IkB = (_PARAM(_a1)*_IkB*_NFkBc-_PARAM(_c6a)*_NFkB_IkB-_PARAM(_a3)*_IKKa*_NFkB_IkB+_PARAM(_e2a)*_NFkB_IkBn) * dt;
  dNFkB_IkBn = (_PARAM(_a1)*KV*_IkBn*_NFkBn-_PARAM(_e2a)*_NFkB_IkBn) * dt;
  dGA20 = (_PARAM(_q1)*_NFkBn*(2-_GA20)-_PARAM(_q2)*_IkBn*_GA20) * dt;
  dGIkB = (_PARAM(_q1)*_NFkBn*(2-_GIkB)-_PARAM(_q2)*_IkBn*_GIkB) * dt;
  dGR = (_PARAM(_q1r)*_NFkBn*(2-_GR)-(_PARAM(_q2rr)+_PARAM(_q2r)*_IkBn)*_GR) * dt;
  dchemt = (_c1rrChemTNF+_c1rChem*_GR-_PARAM(_c3rChem)*_chemt) * dt;
  dchem = (_PARAM(_c4Chem)*_chemt-_PARAM(_c5Chem)*_chem-_PARAM(_e3Chem)*_chem) * dt;
  dTNFt = (_c1rrChemTNF+_c1rTNF*_GR-_PARAM(_c3rTNF)*_TNFt) * dt;
  dTNF = (_PARAM(_c4TNF)*_TNFt-_PARAM(_c5TNF)*_TNF-_PARAM(_e3TNF)*_TNF) * dt;
  dACTt = (_PARAM(_c1rrACT)+_PARAM(_c1r)*_GR-_PARAM(_c3rACT)*_ACTt) * dt;
  dACT = (_PARAM(_c4ACT)*_ACTt-_PARAM(_c5ACT)*_ACT) * dt;
  dIAPt = (_PARAM(_c1rrIAP)+_PARAM(_c1r)*_GR-_PARAM(_c3rIAP)*_IAPt) * dt;
  dIAP = (_PARAM(_c4IAP)*_IAPt-_PARAM(_c5IAP)*_IAP) * dt;

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
  grid.incCCL2(_pos, (_PARAM(_e3Chem) * _chem * dt));
  grid.incCCL5(_pos, (_PARAM(_e3Chem) * _chem * dt));
  grid.incCXCL9(_pos, (2 * _PARAM(_e3Chem) * _chem * dt));

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
  _normalizedACT = _ACT*_PARAM(_c3rACT);
  _IAPt += dIAPt;
  _IAP += dIAP;
  _normalizedIAP = _IAP*_PARAM(_c3rIAP);

  if (_mTNF < 0 || _surfTNFR1 < 0 || _surfBoundTNFR1 < 0 || _surfTNFR2 < 0 || _surfBoundTNFR2 < 0)
    std::cout << "Error: Negative Value of Species in TNF/TNFR dynamics" << std::endl;
  if (_IKKKa < 0 || _IKKn < 0 || _IKKa < 0 || _IKKa < 0 || _IKKi < 0 || _IkBp < 0 || _NFkB_IkBp < 0 ||
      _NFkBc < 0 || _NFkBn < 0 || _A20 < 0 || _A20t < 0 || _IkB < 0 || _IkBn < 0 || _IkBt < 0 ||
      _NFkB_IkB < 0 || _NFkB_IkBn < 0 || _GA20 < 0 || _GIkB < 0 || _GR < 0 || _chemt < 0 || _chem < 0 ||
      _TNFt < 0 || _TNF < 0 || _ACTt < 0 || _ACT < 0 || _IAPt < 0 || _IAP < 0)
    std::cout << "Error: Negative Value of Species in NFkB dynamics" << std::endl;

  //cout << "Debug: Running TNF and NFkB dynamics" << std::endl;
}

void Agent::derivativeTNFandIL10andNFKB(const valarray<double>& vecread, valarray<double>& vecwrite, double dt, GrGrid& grid)
{
  assert(vecread.size() == 41); // Make sure the valarray length is set correctly

  double koff1 = _PARAM(_kOn1) * _PARAM(_KD1);
  double koff2 = _PARAM(_kOn2) * _PARAM(_KD2);
  double Ikoff = _PARAM(_IkOn) * _PARAM(_IkD);
  double IkmRNA;

  // solving for TNF parameters that depend on IL10

  if (_kmRNA > 0)
    {
      IkmRNA = _kmRNA * ((_kSynth/_kmRNA) + ((1.0 - (_kSynth/_kmRNA))/(1.0 + exp((vecread[12] - _PARAM(_LinkRNAGamma))/_PARAM(_LinkRNADelta)))));
    }
  else
    {
      IkmRNA = 0.0;
    }

  // TNF Ordinary Differential Equations

  vecwrite[0] = 0 * dt; // CURRENTLY NOT LINKED WITH IL10 SINCE NFKB DYNAMICS HAVE NOT BEEN LOOKED AT
  vecwrite[1] = (_PARAM(_e3TNF)*_TNF - _kTACE * vecread[1]) * dt;
  vecwrite[2] = (_vTNFR1 - _PARAM(_kOn1) * vecread[8] * vecread[2] + (koff1+KDEG) * vecread[4] - _PARAM(_kT1) * vecread[2] + _PARAM(_kRec1) * vecread[6]) * dt;
  vecwrite[3] = (_vTNFR2 - _PARAM(_kOn2) * vecread[8] * vecread[3] + (koff2+KDEG) * vecread[5] - _PARAM(_kT2) * vecread[3] + _PARAM(_kRec2) * vecread[7]) * dt;
  vecwrite[4] = (_PARAM(_kOn1) * vecread[8] * vecread[2] - (koff1+KDEG) * vecread[4] - _PARAM(_kInt1) * vecread[4]) * dt;
  vecwrite[5] = (_PARAM(_kOn2) * vecread[8] * vecread[3] - (koff2+KDEG) * vecread[5] - _PARAM(_kInt2) * vecread[5] - _PARAM(_kShed) * vecread[5]) * dt;
  vecwrite[6] = (_PARAM(_kInt1) * vecread[4] - _PARAM(_kDeg1) * vecread[6] - _PARAM(_kRec1) * vecread[6]) * dt;
  vecwrite[7] = (_PARAM(_kInt2) * vecread[5] - _PARAM(_kDeg2) * vecread[7] - _PARAM(_kRec2) * vecread[7]) * dt;
  vecwrite[8] = ((DENSITY/NAV) * (_kTACE * vecread[1] - _PARAM(_kOn1) * vecread[8] * vecread[2] + koff1 * vecread[4] - _PARAM(_kOn2) * vecread[8] * vecread[3] + koff2 * vecread[5])) * dt;
  vecwrite[9] = ((DENSITY/NAV) * _PARAM(_kShed) * vecread[5]) * dt;
  // NF-kB dynamics model equations
  vecwrite[10] = (_PARAM(_ka)*vecread[4]*(_PARAM(_KN)-vecread[10])*_PARAM(_kA20)/(_PARAM(_kA20)+vecread[18])-_PARAM(_ki)*vecread[10]) * dt;
  vecwrite[11] = (_PARAM(_k4)*(_PARAM(_KNN)-vecread[11]-vecread[12]-vecread[13])-_PARAM(_k1)*(vecread[10]*vecread[10])*vecread[11]) * dt;
  vecwrite[12] = (_PARAM(_k1)*(vecread[10]*vecread[10])*vecread[11]-_PARAM(_k3)*vecread[12]*(_PARAM(_k2)+vecread[18])/_PARAM(_k2)) * dt;
  vecwrite[13] = (_PARAM(_k3)*vecread[12]*(_PARAM(_k2)+vecread[18])/_PARAM(_k2)-_PARAM(_k4)*vecread[13]) * dt;
  vecwrite[14] = (_PARAM(_a2)*vecread[12]*vecread[20]-_PARAM(_tp)*vecread[14]) * dt;
  vecwrite[15] = (_PARAM(_a3)*vecread[12]*vecread[23]-_PARAM(_tp)*vecread[15]) * dt;
  vecwrite[16] = (_PARAM(_c6a)*vecread[23]-_PARAM(_a1)*vecread[16]*vecread[20]+_PARAM(_tp)*vecread[15]-_PARAM(_i1)*vecread[16]) * dt;
  vecwrite[17] = (_PARAM(_i1)*vecread[16]-_PARAM(_a1)*KV*vecread[21]*vecread[17]) * dt;
  vecwrite[18] = (_PARAM(_c4)*vecread[19]-_PARAM(_c5)*vecread[18]) * dt;
  vecwrite[19] = (_PARAM(_c1)*vecread[25]-_PARAM(_c3)*vecread[19]) * dt;
  vecwrite[20] = (-_PARAM(_a2)*vecread[12]*vecread[20]-_PARAM(_a1)*vecread[20]*vecread[16]+_PARAM(_c4)*vecread[22]-_PARAM(_c5a)*vecread[20]-_PARAM(_i1a)*vecread[20]+_PARAM(_e1a)*vecread[21]) * dt;
  vecwrite[21] = (-_PARAM(_a1)*KV*vecread[21]*vecread[17]+_PARAM(_i1a)*vecread[20]-_PARAM(_e1a)*vecread[21]) * dt;
  vecwrite[22] = (_PARAM(_c1)*vecread[26]-_PARAM(_c3)*vecread[22]) * dt;
  vecwrite[23] = (_PARAM(_a1)*vecread[20]*vecread[16]-_PARAM(_c6a)*vecread[23]-_PARAM(_a3)*vecread[12]*vecread[23]+_PARAM(_e2a)*vecread[24]) * dt;
  vecwrite[24] = (_PARAM(_a1)*KV*vecread[21]*vecread[17]-_PARAM(_e2a)*vecread[24]) * dt;
  vecwrite[25] = (_PARAM(_q1)*vecread[17]*(2-vecread[25])-_PARAM(_q2)*vecread[21]*vecread[25]) * dt;
  vecwrite[26] = (_PARAM(_q1)*vecread[17]*(2-vecread[26])-_PARAM(_q2)*vecread[21]*vecread[26]) * dt;
  vecwrite[27] = (_PARAM(_q1r)*vecread[17]*(2-vecread[27])-(_PARAM(_q2rr)+_PARAM(_q2r)*vecread[21])*vecread[27]) * dt;
  vecwrite[28] = (_c1rrChemTNF+_c1rChem*vecread[27]-_PARAM(_c3rChem)*vecread[28]) * dt;
  vecwrite[29] = (_PARAM(_c4Chem)*vecread[28]-_PARAM(_c5Chem)*vecread[29]-_PARAM(_e3Chem)*vecread[29]) * dt;
  vecwrite[30] = (_c1rrChemTNF+_c1rTNF*vecread[27]-_PARAM(_c3rTNF)*vecread[30]) * dt;
  vecwrite[31] = (_PARAM(_c4TNF)*vecread[30]-_PARAM(_c5TNF)*vecread[31]-_PARAM(_e3TNF)*vecread[31]) * dt;
  vecwrite[32] = (_PARAM(_c1rrACT)+_PARAM(_c1r)*vecread[27]-_PARAM(_c3rACT)*vecread[32]) * dt;
  vecwrite[33] = (_PARAM(_c4ACT)*vecread[32]-_PARAM(_c5ACT)*vecread[33]) * dt;
  vecwrite[34] = (_PARAM(_c1rrIAP)+_PARAM(_c1r)*vecread[27]-_PARAM(_c3rIAP)*vecread[34]) * dt;
  vecwrite[35] = (_PARAM(_c4IAP)*vecread[34]-_PARAM(_c5IAP)*vecread[35]) * dt;
  vecwrite[36] = vecread[33]*_PARAM(_c3rACT);
  vecwrite[37] = vecread[35]*_PARAM(_c3rIAP);

  // IL10 Ordinary Differential Equations
  // sIL10
  vecwrite[38] = (((DENSITY/NAV) * _kISynth) + ((DENSITY/NAV) * (Ikoff * vecread[40] - _PARAM(_IkOn) * vecread[39] * vecread[38]))) * dt;
  // surfIL10R
  vecwrite[39] = (_vIL10R - _PARAM(_IkOn) * vecread[39] * vecread[38] + Ikoff * vecread[40] - _PARAM(_IkT) * vecread[39]) * dt;
  // surfBoundIL10R
  vecwrite[40] = (_PARAM(_IkOn) * vecread[39] * vecread[38] - Ikoff * vecread[40] - _PARAM(_IkInt) * vecread[40]) * dt;

}

void Agent::solveTNFandIL10andNFkB(GrGrid& grid, double dt)
{
  double koff1 = _PARAM(_kOn1) * _PARAM(_KD1);
  double koff2 = _PARAM(_kOn2) * _PARAM(_KD2);
  double Ikoff = _PARAM(_IkOn) * _PARAM(_IkD);

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

  if (_kmRNA > 0)
    {
      IkmRNA = _kmRNA * ((_kSynth/_kmRNA)+ ((1.0 - (_kSynth/_kmRNA))/(1.0 + exp((_surfBoundIL10R - _PARAM(_LinkRNAGamma))/_PARAM(_LinkRNADelta)))));
    }
  else
    {
      IkmRNA = 0.0;
    }

  // end of TNF and IL10 linking


  // CURRENTLY NOT LINKED WITH IL10 SINCE NFKB DYNAMICS HAVE NOT BEEN LOOKED AT

  dmTNF = (_PARAM(_e3TNF)*_TNF - _kTACE * _mTNF) * dt;
  dsurfTNFR1 = (_vTNFR1 - _PARAM(_kOn1) * tnf * _surfTNFR1 + (koff1+KDEG) * _surfBoundTNFR1 - _PARAM(_kT1) * _surfTNFR1 + _PARAM(_kRec1) * _intBoundTNFR1) * dt;
  dsurfTNFR2 = (_vTNFR2 - _PARAM(_kOn2) * tnf * _surfTNFR2 + (koff2+KDEG) * _surfBoundTNFR2 - _PARAM(_kT2) * _surfTNFR2 + _PARAM(_kRec2) * _intBoundTNFR2) * dt;
  dsurfBoundTNFR1 = (_PARAM(_kOn1) * tnf * _surfTNFR1 - (koff1+KDEG) * _surfBoundTNFR1 - _PARAM(_kInt1) * _surfBoundTNFR1) * dt;
  dsurfBoundTNFR2 = (_PARAM(_kOn2) * tnf * _surfTNFR2 - (koff2+KDEG) * _surfBoundTNFR2 - _PARAM(_kInt2) * _surfBoundTNFR2 - _PARAM(_kShed) * _surfBoundTNFR2) * dt;
  dintBoundTNFR1 = (_PARAM(_kInt1) * _surfBoundTNFR1 - _PARAM(_kDeg1) * _intBoundTNFR1 - _PARAM(_kRec1) * _intBoundTNFR1) * dt;
  dintBoundTNFR2 = (_PARAM(_kInt2) * _surfBoundTNFR2 - _PARAM(_kDeg2) * _intBoundTNFR2 - _PARAM(_kRec2) * _intBoundTNFR2) * dt;
  dsTNF = ((DENSITY/NAV) * (_kTACE * _mTNF - _PARAM(_kOn1) * tnf * _surfTNFR1 + koff1 * _surfBoundTNFR1 - _PARAM(_kOn2) * tnf * _surfTNFR2 + koff2 * _surfBoundTNFR2)) * dt;
  dshedTNFR2 = ((DENSITY/NAV) * _PARAM(_kShed) * _surfBoundTNFR2) * dt;

  // NF-kB dynamics model equations
  dIKKKa = (_PARAM(_ka)*_surfBoundTNFR1*(_PARAM(_KN)-_IKKKa)*_PARAM(_kA20)/(_PARAM(_kA20)+_A20)-_PARAM(_ki)*_IKKKa) * dt;
  dIKKn = (_PARAM(_k4)*(_PARAM(_KNN)-_IKKn-_IKKa-_IKKi)-_PARAM(_k1)*(_IKKKa*_IKKKa)*_IKKn) * dt;
  dIKKa = (_PARAM(_k1)*(_IKKKa*_IKKKa)*_IKKn-_PARAM(_k3)*_IKKa*(_PARAM(_k2)+_A20)/_PARAM(_k2)) * dt;
  dIKKi = (_PARAM(_k3)*_IKKa*(_PARAM(_k2)+_A20)/_PARAM(_k2)-_PARAM(_k4)*_IKKi) * dt;
  dIkBp = (_PARAM(_a2)*_IKKa*_IkB-_PARAM(_tp)*_IkBp) * dt;
  dNFkB_IkBp = (_PARAM(_a3)*_IKKa*_NFkB_IkB-_PARAM(_tp)*_NFkB_IkBp) * dt;
  dNFkBc = (_PARAM(_c6a)*_NFkB_IkB-_PARAM(_a1)*_NFkBc*_IkB+_PARAM(_tp)*_NFkB_IkBp-_PARAM(_i1)*_NFkBc) * dt;
  dNFkBn = (_PARAM(_i1)*_NFkBc-_PARAM(_a1)*KV*_IkBn*_NFkBn) * dt;
  dA20 = (_PARAM(_c4)*_A20t-_PARAM(_c5)*_A20) * dt;
  dA20t = (_PARAM(_c1)*_GA20-_PARAM(_c3)*_A20t) * dt;
  dIkB = (-_PARAM(_a2)*_IKKa*_IkB-_PARAM(_a1)*_IkB*_NFkBc+_PARAM(_c4)*_IkBt-_PARAM(_c5a)*_IkB-_PARAM(_i1a)*_IkB+_PARAM(_e1a)*_IkBn) * dt;
  dIkBn = (-_PARAM(_a1)*KV*_IkBn*_NFkBn+_PARAM(_i1a)*_IkB-_PARAM(_e1a)*_IkBn) * dt;
  dIkBt = (_PARAM(_c1)*_GIkB-_PARAM(_c3)*_IkBt) * dt;
  dNFkB_IkB = (_PARAM(_a1)*_IkB*_NFkBc-_PARAM(_c6a)*_NFkB_IkB-_PARAM(_a3)*_IKKa*_NFkB_IkB+_PARAM(_e2a)*_NFkB_IkBn) * dt;
  dNFkB_IkBn = (_PARAM(_a1)*KV*_IkBn*_NFkBn-_PARAM(_e2a)*_NFkB_IkBn) * dt;
  dGA20 = (_PARAM(_q1)*_NFkBn*(2-_GA20)-_PARAM(_q2)*_IkBn*_GA20) * dt;
  dGIkB = (_PARAM(_q1)*_NFkBn*(2-_GIkB)-_PARAM(_q2)*_IkBn*_GIkB) * dt;
  dGR = (_PARAM(_q1r)*_NFkBn*(2-_GR)-(_PARAM(_q2rr)+_PARAM(_q2r)*_IkBn)*_GR) * dt;
  dchemt = (_c1rrChemTNF+_c1rChem*_GR-_PARAM(_c3rChem)*_chemt) * dt;
  dchem = (_PARAM(_c4Chem)*_chemt-_PARAM(_c5Chem)*_chem-_PARAM(_e3Chem)*_chem) * dt;
  dTNFt = (_c1rrChemTNF+_c1rTNF*_GR-_PARAM(_c3rTNF)*_TNFt) * dt;
  dTNF = (_PARAM(_c4TNF)*_TNFt-_PARAM(_c5TNF)*_TNF-_PARAM(_e3TNF)*_TNF) * dt;
  dACTt = (_PARAM(_c1rrACT)+_PARAM(_c1r)*_GR-_PARAM(_c3rACT)*_ACTt) * dt;
  dACT = (_PARAM(_c4ACT)*_ACTt-_PARAM(_c5ACT)*_ACT) * dt;
  dIAPt = (_PARAM(_c1rrIAP)+_PARAM(_c1r)*_GR-_PARAM(_c3rIAP)*_IAPt) * dt;
  dIAP = (_PARAM(_c4IAP)*_IAPt-_PARAM(_c5IAP)*_IAP) * dt;

  // IL10 differential equations
  dsIL10 = ((DENSITY/NAV) * _kISynth + ((DENSITY/NAV) * (Ikoff * _surfBoundIL10R - _PARAM(_IkOn) * _surfIL10R * il10))) * dt;
  dsurfIL10R = (_vIL10R - _PARAM(_IkOn) * _surfIL10R * il10 + Ikoff * _surfBoundIL10R - _PARAM(_IkT) * _surfIL10R) * dt;
  dsurfBoundIL10R = (_PARAM(_IkOn) * _surfIL10R * il10 - Ikoff * _surfBoundIL10R - _PARAM(_IkInt) * _surfBoundIL10R) * dt;
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
  grid.incCCL2(_pos, (_PARAM(_e3Chem) * _chem * dt));
  grid.incCCL5(_pos, (_PARAM(_e3Chem) * _chem * dt));
  grid.incCXCL9(_pos,  (2 * _PARAM(_e3Chem) * _chem * dt));

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
  _normalizedACT = _ACT*_PARAM(_c3rACT);
  _IAPt += dIAPt;
  _IAP += dIAP;
  _normalizedIAP = _IAP*_PARAM(_c3rIAP);

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

#endif

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
  //dIKKKa = (_PARAM(_ka)*_surfBoundTNFR1*(_PARAM(_KN)-_IKKKa)*_PARAM(_kA20)/(_PARAM(_kA20)+_A20)-_PARAM(_ki)*_IKKKa) * dt;
  dIKKKa = 0.0;
  //dIKKn = (_PARAM(_k4)*(_PARAM(_KNN)-_IKKn-_IKKa-_IKKi)-_PARAM(_k1)*(_IKKKa*_IKKKa)*_IKKn) * dt;
  dIKKn = 0.0;
  //dIKKa = (_PARAM(_k1)*(_IKKKa*_IKKKa)*_IKKn-_PARAM(_k3)*_IKKa*(_PARAM(_k2)+_A20)/_PARAM(_k2)) * dt;
  dIKKa = 0.0;
  //dIKKi = (_PARAM(_k3)*_IKKa*(_PARAM(_k2)+_A20)/_PARAM(_k2)-_PARAM(_k4)*_IKKi) * dt;
  dIKKi = 0.0;
  //dIkBp = (_PARAM(_a2)*_IKKa*_IkB-_PARAM(_tp)*_IkBp) * dt;
  dIkBp = 0.0;
  //dNFkB_IkBp = (_PARAM(_a3)*_IKKa*_NFkB_IkB-_PARAM(_tp)*_NFkB_IkBp) * dt;
  dNFkB_IkBp = 0.0;
  dNFkBc = (_PARAM(_c6a)*_NFkB_IkB-_PARAM(_a1)*_NFkBc*_IkB+_PARAM(_tp)*_NFkB_IkBp-_PARAM(_i1)*_NFkBc) * dt;
  dNFkBn = (_PARAM(_i1)*_NFkBc-_PARAM(_a1)*KV*_IkBn*_NFkBn) * dt;
  dA20 = (_PARAM(_c4)*_A20t-_PARAM(_c5)*_A20) * dt;
  dA20t = (_PARAM(_c1)*_GA20-_PARAM(_c3)*_A20t) * dt;
  dIkB = (-_PARAM(_a2)*_IKKa*_IkB-_PARAM(_a1)*_IkB*_NFkBc+_PARAM(_c4)*_IkBt-_PARAM(_c5a)*_IkB-_PARAM(_i1a)*_IkB+_PARAM(_e1a)*_IkBn) * dt;
  dIkBn = (-_PARAM(_a1)*KV*_IkBn*_NFkBn+_PARAM(_i1a)*_IkB-_PARAM(_e1a)*_IkBn) * dt;
  dIkBt = (_PARAM(_c1)*_GIkB-_PARAM(_c3)*_IkBt) * dt;
  dNFkB_IkB = (_PARAM(_a1)*_IkB*_NFkBc-_PARAM(_c6a)*_NFkB_IkB-_PARAM(_a3)*_IKKa*_NFkB_IkB+_PARAM(_e2a)*_NFkB_IkBn) * dt;
  dNFkB_IkBn = (_PARAM(_a1)*KV*_IkBn*_NFkBn-_PARAM(_e2a)*_NFkB_IkBn) * dt;
  dGA20 = (_PARAM(_q1)*_NFkBn*(2-_GA20)-_PARAM(_q2)*_IkBn*_GA20) * dt;
  dGIkB = (_PARAM(_q1)*_NFkBn*(2-_GIkB)-_PARAM(_q2)*_IkBn*_GIkB) * dt;
  dGR = (_PARAM(_q1r)*_NFkBn*(2-_GR)-(_PARAM(_q2rr)+_PARAM(_q2r)*_IkBn)*_GR) * dt;
  //dchemt = (_c1rrChemTNF+_c1rChem*_GR-_PARAM(_c3rChem)*_chemt) * dt;
  dchemt = 0.0;
  //dchem = (_PARAM(_c4Chem)*_chemt-_PARAM(_c5Chem)*_chem-_PARAM(_e3Chem)*_chem) * dt;
  dchem = 0.0;
  //dTNFt = (_c1rrChemTNF+_c1rTNF*_GR-_PARAM(_c3rTNF)*_TNFt) * dt;
  dTNFt = 0.0;
  //dTNF = (_PARAM(_c4TNF)*_TNFt-_PARAM(_c5TNF)*_TNF-_PARAM(_e3TNF)*_TNF) * dt;
  dTNF = 0.0;
  //dACTt = (_PARAM(_c1rrACT)+_PARAM(_c1r)*_GR-_PARAM(_c3rACT)*_ACTt) * dt;
  dACTt = 0.0;
  //dACT = (_PARAM(_c4ACT)*_ACTt-_PARAM(_c5ACT)*_ACT) * dt;
  dACT = 0.0;
  //dIAPt = (_PARAM(_c1rrIAP)+_PARAM(_c1r)*_GR-_PARAM(_c3rIAP)*_IAPt) * dt;
  dIAPt = 0.0;
  //dIAP = (_PARAM(_c4IAP)*_IAPt-_PARAM(_c5IAP)*_IAP) * dt;
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
  _normalizedACT = _ACT*_PARAM(_c3rACT);
  _IAPt += dIAPt;
  _IAP += dIAP;
  _normalizedIAP = _IAP*_PARAM(_c3rIAP);
  /*
   if (_IKKKa < 0 || _IKKn < 0 || _IKKa < 0 || _IKKa < 0 || _IKKi < 0 || _IkBp < 0 || _NFkB_IkBp < 0 ||
   _NFkBc < 0 || _NFkBn < 0 || _A20 < 0 || _A20t < 0 || _IkB < 0 || _IkBn < 0 || _IkBt < 0 ||
   _NFkB_IkB < 0 || _NFkB_IkBn < 0 || _GA20 < 0 || _GIkB < 0 || _GR < 0 || _chemt < 0 || _chem < 0 ||
   _TNFt < 0 || _TNF < 0 || _ACTt < 0 || _ACT < 0 || _IAPt < 0 || _IAP < 0)
   std::cout << "Error: Negative Value of Species in NFkB dynamics" << std::endl;
   */
}



bool Agent::TNFinducedApoptosis(GrGrid &grid, bool tnfrDynamics, bool nfkbDynamics)
{
  bool Threshold, Probability;

  if (!nfkbDynamics && tnfrDynamics)
    {

      Scalar testintBoundTNFR1 = std::min(_intBoundTNFR1, _PARAM(_saturationApoptosisTNF_Molecular));

      Threshold = intCompareGT(_intBoundTNFR1, _PARAM(_thresholdApoptosisTNF_Molecular));
      Probability = intCompareGT(1 - exp(-_PARAM(_kApoptosis_Molecular) * (testintBoundTNFR1 - _PARAM(_thresholdApoptosisTNF_Molecular))), g_Rand.getReal());

    }
  else if (nfkbDynamics)
    {
      double nfkb_adjusted_k_apoptosis = _PARAM(_kApoptosis_NFkB_Molecular) * (_PARAM(_kIAP)/(_PARAM(_kIAP) + _normalizedIAP));

      Threshold = intCompareGT(_intBoundTNFR1, _PARAM(_thresholdApoptosisTNF_Molecular));
      Probability = intCompareGT(1 - exp(-nfkb_adjusted_k_apoptosis * (_intBoundTNFR1 - _PARAM(_thresholdApoptosisTNF_Molecular))), g_Rand.getReal());

    }
  else if (!tnfrDynamics)
    {
      double tnfBoundFraction = grid.TNF(_pos) / (grid.TNF(_pos) + _PARAM(_KD1) * 48.16e11);

      Threshold = intCompareGT(tnfBoundFraction, _PARAM(_thresholdApoptosisTNF));
      Probability = intCompareGT(1 - exp(-_PARAM(_kApoptosis) * (tnfBoundFraction - _PARAM(_thresholdApoptosisTNF))), g_Rand.getReal());

    }
  else
    {
      throw std::runtime_error("Invalid TNF Induced Apoptosis Calculation -- Criteria Not Met For Molecular/Non-Molecular Options");
    }

  return Threshold && Probability;

}

bool Agent::TNFinducedNFkB(GrGrid &grid, bool tnfrDynamics, bool nfkbDynamics)
{
  bool Threshold, Probability;

  if (!nfkbDynamics && tnfrDynamics)
    {
      Scalar testsurfBoundTNFR1;

      testsurfBoundTNFR1 = std::min(_surfBoundTNFR1, (_PARAM(Mac_saturationFracNFkBTNF_Molecular) * (_surfTNFR1 + _surfBoundTNFR1)));

      Threshold = intCompareGT(_surfBoundTNFR1, _PARAM(Mac_thresholdNFkBTNF_Molecular));
      Probability = intCompareGT(1 - exp(-_PARAM(Mac_kNFkB_Molecular) * (testsurfBoundTNFR1 - _PARAM(Mac_thresholdNFkBTNF_Molecular))), g_Rand.getReal());
    }


  else if (nfkbDynamics)
    {
      Threshold = intCompareGT(_normalizedACT, _PARAM(_actThreshold));
      Probability = intCompareGT(1 - exp(-_PARAM(_activationRate) * (_normalizedACT - _PARAM(_actThreshold))), g_Rand.getReal());
    }


  else if (!tnfrDynamics)
    {
      double tnfBoundFraction = grid.TNF(_pos) / (grid.TNF(_pos) + _PARAM(_KD1) * 48.16e11);

      Threshold = intCompareGT(tnfBoundFraction, _PARAM(Mac_thresholdNFkBTNF));
      Probability = intCompareGT(1 - exp(-_PARAM(Mac_kNFkB) * (tnfBoundFraction - _PARAM(Mac_thresholdNFkBTNF))), g_Rand.getReal());
    }

  else
    {
      throw std::runtime_error("Invalid TNF Induced NFkB Calculation -- Criteria Not Met For Molecular/Non-Molecular Options");
    }

  return Threshold && Probability;

}


void Agent::calcIkmRNA(GrGrid& grid, double& kmRNA, double ksynth, bool il10rDynamics)
{
  double eqsurfBoundIL10R;

  if(il10rDynamics)
    eqsurfBoundIL10R = _surfBoundIL10R;
  else
    eqsurfBoundIL10R = ((grid.il10(_pos)/(NAV * VOL)) * _surfIL10R) / (_PARAM(_IkD) + (grid.il10(_pos)/(NAV * VOL)));

  kmRNA = (kmRNA) * ((ksynth/(kmRNA)) + ((1.0 - (ksynth/(kmRNA)))/(1.0 + exp(((eqsurfBoundIL10R - _PARAM(_LinkRNAGamma))/_PARAM(_LinkRNADelta))))));

}

void Agent::solveDegradation(GrGrid& grid, double dt, bool tnfrDynamics, bool il10rDynamics, Scalar meanTNFR1, Scalar iIL10R)
{
  if (!tnfrDynamics)
    {

      // simulate the effect of TNF internalization by cells in the form of degradation. Only for TNF
      double dtnf;
      double tnf = grid.TNF(_pos);
      dtnf = -_PARAM(_kInt1) * (tnf / (tnf + _PARAM(_KD1) * NAV * VOL)) * meanTNFR1 * dt * 0.4;
      grid.incTNF(_pos, dtnf);
    }

  if (!il10rDynamics)
    {

      double dil10;
      double il10 = grid.il10(_pos);

      // simulate the effect of IL10 internalization in the form of degradation. Only for IL10
      dil10 = -_PARAM(_IkInt) * (il10 / (il10 + _PARAM(_IkD) * NAV * VOL)) * iIL10R * dt * _PARAM(_Imod);
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

  double min = _PARAM(_minChemotaxis);
  double max = _PARAM(_maxChemotaxis);

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
      for (k = 0; k < 9 && cumProb[k] < r; k++) ;
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

#if 0
void Agent::serialize(std::ostream& out) const
{
  assert(out.good());

  Serialization::writeHeader(out, Agent::_ClassName);

#define P(type, name, ival, desc) \
    out << _##name << std::endl;
  AGENT_PROPS
#undef P
  out << _lasttimestep << std::endl;

  out << _initvector.size() << endl;

  for (size_t jj = 0; jj < _initvector.size(); jj++)
    {
      out << _initvector[jj] << std::endl;
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

#define P(type, name, ival, desc) \
    in >> _##name;
  AGENT_PROPS
#undef P
  in >> _lasttimestep;

  size_t vectorSize;
  in >> vectorSize;

  _initvector.resize(vectorSize, 0.0);


  for (size_t kk = 0; kk < vectorSize; kk++)
    {
      in >> _initvector[kk];
    }

  if (!Serialization::readFooter(in, Agent::_ClassName))
    {
      exit(1);
    }
}
#endif

/*virtual*/ void LungFunc::operator()(const ODESolvers::ODEState& vecread, double t, ODESolvers::Derivative& vecwrite, void* params) const
{
  const Agent* const agent = ((Params_t*)params)->agent;  //params will always be passed a Params_t, cast it and use it as such
  const GrGrid& grid = *(((Params_t*)params)->grid);
  const double koff1 = _PARAM(_kOn1) * _PARAM(_KD1);
  const double koff2 = _PARAM(_kOn2) * _PARAM(_KD2);
  double il10 = grid.il10(agent->getPosition()) / (NAV * VOL);

  if(_PARAM(_TNFdynamics))
    {
      // solving for TNF parameters that depend on IL10
      double eqsurfBoundIL10R;
      if(_PARAM(_IL10dynamics))
        eqsurfBoundIL10R = vecread[il10offset+2];
      else
        eqsurfBoundIL10R = (il10 * agent->getmeanIL10R()) / (_PARAM(_IkD) + il10);

      double IkmRNA;
      if (agent->getkmRNA() > 0)
        IkmRNA = agent->getkmRNA() * ((agent->getkSynth()/agent->getkmRNA()) + ((1.0 - (agent->getkSynth()/agent->getkmRNA()))/(1.0 + exp((eqsurfBoundIL10R - _PARAM(_LinkRNAGamma))/_PARAM(_LinkRNADelta)))));
      else
        IkmRNA = 0.0;

//    double IkmRNA = agent->getkmRNA();

      // TNF Ordinary Differential Equations
      // mTNFRNA
      vecwrite[0] = (IkmRNA - _PARAM(_kTrans) * vecread[0]);
      // mTNF
      vecwrite[1] = (_PARAM(_kTrans) * vecread[0] - agent->getkTACE() * vecread[1]);
      // surfTNFR1
      vecwrite[2] = (agent->getvTNFR1() - _PARAM(_kOn1) * (vecread[8] * 1.0/(NAV * VOL)) * vecread[2] + koff1 * vecread[4] - _PARAM(_kT1) * vecread[2] + _PARAM(_kRec1) * vecread[6]);
      // surfTNFR2
      vecwrite[3] = (agent->getvTNFR2() - _PARAM(_kOn2) * (vecread[8] * 1.0/(NAV * VOL)) * vecread[3] + koff2 * vecread[5] - _PARAM(_kT2) * vecread[3] + _PARAM(_kRec2) * vecread[7]);
      // surfBoundTNFR1
      vecwrite[4] = (_PARAM(_kOn1) * (vecread[8] * 1.0/(NAV * VOL)) * vecread[2] - koff1 * vecread[4] - _PARAM(_kInt1) * vecread[4]);
      // surfBoundTNFR2
      vecwrite[5] = (_PARAM(_kOn2) * (vecread[8] * 1.0/(NAV * VOL)) * vecread[3] - koff2 * vecread[5] - _PARAM(_kInt2) * vecread[5] - _PARAM(_kShed) * vecread[5]);
      // intBoundTNFR1
      vecwrite[6] = (_PARAM(_kInt1) * vecread[4] - _PARAM(_kDeg1) * vecread[6] - _PARAM(_kRec1) * vecread[6]);
      // intBoundTNFR2
      vecwrite[7] = (_PARAM(_kInt2) * vecread[5] - _PARAM(_kDeg2) * vecread[7] - _PARAM(_kRec2) * vecread[7]);
      // sTNF
      vecwrite[8] = (((DENSITY/NAV) * (agent->getkTACE() * vecread[1] - _PARAM(_kOn1) * (vecread[8] * 1.0/(NAV * VOL)) * vecread[2] + koff1 * vecread[4] - _PARAM(_kOn2) * (vecread[8] * 1.0/(NAV * VOL)) * vecread[3] + koff2 * vecread[5]))) * (NAV * VOL);
      // shedTNFR2
      vecwrite[9] = (((DENSITY/NAV) * _PARAM(_kShed) * vecread[5])) * (NAV * VOL);
    }
  if(_PARAM(_IL10dynamics))
    il10deriv(vecread, t, vecwrite, ((Params_t*)params));
}

inline void LungFunc::il10deriv(const ODESolvers::ODEState& vecread, double /*t*/, ODESolvers::Derivative& vecwrite, Params_t* params) const
{
  const double Ikoff = _PARAM(_IkOn) * _PARAM(_IkD);
  // IL10 Ordinary Differential Equations
  // sIL10
  vecwrite[il10offset+0] = ((((DENSITY/NAV) * params->agent->getkISynth()) + ((DENSITY/NAV) * (Ikoff * vecread[il10offset+2] - _PARAM(_IkOn) * vecread[il10offset+1] * (vecread[il10offset+0] * 1.0/(NAV * VOL)))))) * (NAV * VOL);
  // surfIL10R
  vecwrite[il10offset+1] = (params->agent->getvIL10R() - _PARAM(_IkOn) * vecread[il10offset+1] * (vecread[il10offset+0] * 1.0/(NAV * VOL)) + Ikoff * vecread[il10offset+2] - _PARAM(_IkT) * vecread[il10offset+1]);
  // surfBoundIL10R
  vecwrite[il10offset+2] = (_PARAM(_IkOn) * vecread[il10offset+1] * (vecread[il10offset+0] * 1.0/(NAV * VOL)) - Ikoff * vecread[il10offset+2] - _PARAM(_IkInt) * vecread[il10offset+2]);
}

void NFKBFunc::operator()(const ODESolvers::ODEState& vecread, double t, ODESolvers::Derivative& vecwrite, void* params) const
{
  const Agent* const agent = ((Params_t*)params)->agent;  //params will always be passed a Params_t, cast it and use it as such
  //const GrGrid& grid = *(((Params_t*)params)->grid);
  const double koff1 = _PARAM(_kOn1) * _PARAM(_KD1);
  const double koff2 = _PARAM(_kOn2) * _PARAM(_KD2);
  //double il10 = grid.il10(agent->getPosition()) / (NAV * VOL);
  assert(_PARAM(_NFkBdynamics) && "NFKB NOT ENABLED!");
  assert(agent->NFkBCapable() && "NFKB ON NON-NFKB CAPABLE AGENT");
  vecwrite[0] = 0; // CURRENTLY NOT LINKED WITH IL10 SINCE NFKB DYNAMICS HAVE NOT BEEN LOOKED AT
  vecwrite[1] = (_PARAM(_e3TNF)*agent->getTNF() - agent->getkTACE() * vecread[1]);
  vecwrite[2] = (agent->getvTNFR1() - _PARAM(_kOn1) * (vecread[8]/(NAV * VOL)) * vecread[2] + (koff1+KDEG) * vecread[4] - _PARAM(_kT1) * vecread[2] + _PARAM(_kRec1) * vecread[6]);
  vecwrite[3] = (agent->getvTNFR2() - _PARAM(_kOn2) * (vecread[8]/(NAV * VOL)) * vecread[3] + (koff2+KDEG) * vecread[5] - _PARAM(_kT2) * vecread[3] + _PARAM(_kRec2) * vecread[7]);
  vecwrite[4] = (_PARAM(_kOn1) * (vecread[8]/(NAV * VOL)) * vecread[2] - (koff1+KDEG) * vecread[4] - _PARAM(_kInt1) * vecread[4]);
  vecwrite[5] = (_PARAM(_kOn2) * (vecread[8]/(NAV * VOL)) * vecread[3] - (koff2+KDEG) * vecread[5] - _PARAM(_kInt2) * vecread[5] - _PARAM(_kShed) * vecread[5]);
  vecwrite[6] = (_PARAM(_kInt1) * vecread[4] - _PARAM(_kDeg1) * vecread[6] - _PARAM(_kRec1) * vecread[6]);
  vecwrite[7] = (_PARAM(_kInt2) * vecread[5] - _PARAM(_kDeg2) * vecread[7] - _PARAM(_kRec2) * vecread[7]);
  vecwrite[8] = (((DENSITY/NAV) * (agent->getkTACE() * vecread[1] - _PARAM(_kOn1) * (vecread[8]/(NAV * VOL)) * vecread[2] + koff1 * vecread[4] - _PARAM(_kOn2) * (vecread[8]/(NAV * VOL)) * vecread[3] + koff2 * vecread[5]))) * (NAV * VOL);
  vecwrite[9] = (((DENSITY/NAV) * _PARAM(_kShed) * vecread[5])) * (NAV * VOL);
  // NF-kB dynamics model equations
  vecwrite[10] = (_PARAM(_ka)*vecread[4]*(_PARAM(_KN)-vecread[10])*_PARAM(_kA20)/(_PARAM(_kA20)+vecread[18])-_PARAM(_ki)*vecread[10]);
  vecwrite[11] = (_PARAM(_k4)*(_PARAM(_KNN)-vecread[11]-vecread[12]-vecread[13])-_PARAM(_k1)*(vecread[10]*vecread[10])*vecread[11]);
  vecwrite[12] = (_PARAM(_k1)*(vecread[10]*vecread[10])*vecread[11]-_PARAM(_k3)*vecread[12]*(_PARAM(_k2)+vecread[18])/_PARAM(_k2));
  vecwrite[13] = (_PARAM(_k3)*vecread[12]*(_PARAM(_k2)+vecread[18])/_PARAM(_k2)-_PARAM(_k4)*vecread[13]);
  vecwrite[14] = (_PARAM(_a2)*vecread[12]*vecread[20]-_PARAM(_tp)*vecread[14]);
  vecwrite[15] = (_PARAM(_a3)*vecread[12]*vecread[23]-_PARAM(_tp)*vecread[15]);
  vecwrite[16] = (_PARAM(_c6a)*vecread[23]-_PARAM(_a1)*vecread[16]*vecread[20]+_PARAM(_tp)*vecread[15]-_PARAM(_i1)*vecread[16]);
  vecwrite[17] = (_PARAM(_i1)*vecread[16]-_PARAM(_a1)*KV*vecread[21]*vecread[17]);
  vecwrite[18] = (_PARAM(_c4)*vecread[19]-_PARAM(_c5)*vecread[18]);
  vecwrite[19] = (_PARAM(_c1)*vecread[25]-_PARAM(_c3)*vecread[19]);
  vecwrite[20] = (-_PARAM(_a2)*vecread[12]*vecread[20]-_PARAM(_a1)*vecread[20]*vecread[16]+_PARAM(_c4)*vecread[22]-_PARAM(_c5a)*vecread[20]-_PARAM(_i1a)*vecread[20]+_PARAM(_e1a)*vecread[21]);
  vecwrite[21] = (-_PARAM(_a1)*KV*vecread[21]*vecread[17]+_PARAM(_i1a)*vecread[20]-_PARAM(_e1a)*vecread[21]);
  vecwrite[22] = (_PARAM(_c1)*vecread[26]-_PARAM(_c3)*vecread[22]);
  vecwrite[23] = (_PARAM(_a1)*vecread[20]*vecread[16]-_PARAM(_c6a)*vecread[23]-_PARAM(_a3)*vecread[12]*vecread[23]+_PARAM(_e2a)*vecread[24]);
  vecwrite[24] = (_PARAM(_a1)*KV*vecread[21]*vecread[17]-_PARAM(_e2a)*vecread[24]);
  vecwrite[25] = (_PARAM(_q1)*vecread[17]*(2-vecread[25])-_PARAM(_q2)*vecread[21]*vecread[25]);
  vecwrite[26] = (_PARAM(_q1)*vecread[17]*(2-vecread[26])-_PARAM(_q2)*vecread[21]*vecread[26]);
  vecwrite[27] = (_PARAM(_q1r)*vecread[17]*(2-vecread[27])-(_PARAM(_q2rr)+_PARAM(_q2r)*vecread[21])*vecread[27]);
  vecwrite[28] = (agent->getc1rrChemTNF()+agent->getc1rChem()*vecread[27]-_PARAM(_c3rChem)*vecread[28]);
  vecwrite[29] = (_PARAM(_c4Chem)*vecread[28]-_PARAM(_c5Chem)*vecread[29]-_PARAM(_e3Chem)*vecread[29]);
  vecwrite[30] = (agent->getc1rrChemTNF()+agent->getc1rTNF()*vecread[27]-_PARAM(_c3rTNF)*vecread[30]);
  vecwrite[31] = (_PARAM(_c4TNF)*vecread[30]-_PARAM(_c5TNF)*vecread[31]-_PARAM(_e3TNF)*vecread[31]);
  vecwrite[32] = (_PARAM(_c1rrACT)+_PARAM(_c1r)*vecread[27]-_PARAM(_c3rACT)*vecread[32]);
  vecwrite[33] = (_PARAM(_c4ACT)*vecread[32]-_PARAM(_c5ACT)*vecread[33]);
  vecwrite[34] = (_PARAM(_c1rrIAP)+_PARAM(_c1r)*vecread[27]-_PARAM(_c3rIAP)*vecread[34]);
  vecwrite[35] = (_PARAM(_c4IAP)*vecread[34]-_PARAM(_c5IAP)*vecread[35]);
#if 0 //Not necassary, done on copy
  vecwrite[36] = vecread[33]*_PARAM(_c3rACT); //These should not be multiplied by dt (stepper will do so automatically)
  vecwrite[37] = vecread[35]*_PARAM(_c3rIAP);
#endif
  if(_PARAM(_IL10dynamics))
    il10deriv(vecread, t, vecwrite, ((Params_t*)params));
}

#if 1 //Possible implementation of adaptive 2 cell ode integration
void Agent::solveMolecularScaleAdaptive(GrGrid& grid, double t, double dt, ODESolvers::ODEMethod method)
{
  Agent* other = grid.agent(getPosition(), 0) == this ? grid.agent(getPosition(), 1) : grid.agent(getPosition(), 0);
  if(other)
    {
      assert(getODEstatus() == other->getODEstatus());
      adaptiveODE2Cell(other, grid, t, dt, method);
      other->setODEstatus(true);
    }
  else
    adaptiveODE(grid, t, dt, method);
  setODEstatus(true);
}

#define SAFETY (0.95)
#define PGROW (-0.2)
#define PSHRINK (-0.25)
#define ERRCON (1.89e-4)
// Define Accuracy
#define ACCURACY (1e-5)
#define TINY (1e-30)
// Define Max # of Steps
#define MAXSTEP (10000)
// Define Min Step Size
#define MIN_STEP_SIZE (0.01)

#define EPS (1e-5)

#define ABSTOL (1e-6)
#define RELTOL (1e-5)
#define MAXSCALE (10.0)
#define MINSCALE (0.2)

void Agent::adaptiveODE(GrGrid& grid, double t, double dt, ODESolvers::ODEMethod method)
{
  writeValarrayFromMembers(grid, _initvector);  //With some pointer magic, should be able to remove this...
  LungFunc& fn = *getDerivFunc();
  const size_t  dim = fn.dim(); // Avoid calling fn.dim() at all costs
  valarray<double> error(dim);
  static LungFunc::Params_t params;
  params.agent = this;
  params.grid = &grid;
  ODESolvers::Stepper* s1 = getStepper(method);
  ODESolvers::AdaptiveStepper stepper(dim, s1, EPS, TINY, MAXSTEP, PSHRINK, PGROW, SAFETY, ABSTOL, RELTOL, MAXSCALE, MINSCALE); // Just use the adaptive stepper from ODESolvers, not very efficient
  stepper.step(_initvector, fn, t, dt, _lasttimestep, error, (void*)&params);
  writeMembersFromValarray(grid, _initvector);
}

void Agent::adaptiveODE2Cell(Agent* other, GrGrid& grid, double t, double dt, ODESolvers::ODEMethod method)
{
  //Rewrite the adaptive stepper method with two sets of odes
  assert(other && other->getPosition() == this->getPosition());
  writeValarrayFromMembers(grid, _initvector);
  other->writeValarrayFromMembers(grid, other->_initvector);
  static LungFunc::Params_t params;
  params.grid = &grid;
  ODESolvers::Stepper* s1 = getStepper(method);
  ODESolvers::Stepper* s2 = other->getStepper(method);
  LungFunc& fn1 = *getDerivFunc();
  LungFunc& fn2 = *(other->getDerivFunc());
  valarray<double> dy1(fn1.dim()), err1(fn1.dim()), yscal1(fn1.dim()), tmp1(_initvector);
  valarray<double> dy2(fn2.dim()), err2(fn2.dim()), yscal2(fn2.dim()), tmp2(other->_initvector);
  const double t1 = t, t2 = t+dt;
  static double h = dt; //Static here to adjust for next set of agents
  //Time management
  for(size_t i=0; i<MAXSTEP; i++)
    {
      {
        //fn(i,t,dy,params); yscal = abs(i) + abs(dy)*dt + TINY;
        params.agent = this;
        fn1(_initvector, t, dy1, (void*)&params);
        yscal1 = abs(tmp1) + abs(dy1)*dt + TINY;
        params.agent = other;
        fn2(other->_initvector, t, dy2, (void*)&params);
        yscal2 = abs(tmp2) + abs(dy2)*dt + TINY;
      }
      if((t+dt-t2)*(t+dt-t1) > 0.0) //Finish up the last timestep
        dt = t2 - t;
      //Step management
      for(;;)
        {
          {
            //stepper->step(tmp, fn, t, h, err, params);
            params.agent = this;
            s1->step(tmp1, fn1, t, h, _lasttimestep, err1, (void*)&params);

            if(_PARAM(_TNFdynamics)) //Copy to the second agent
              tmp2[fn2.tnfidx()] = tmp1[fn1.tnfidx()];
            if(_PARAM(_IL10dynamics))
              tmp2[fn2.il10idx()] = tmp1[fn1.il10idx()];

            params.agent = other;
            s2->step(tmp2, fn2, t, h, _lasttimestep, err2, (void*)&params);

            if(_PARAM(_TNFdynamics)) //Copy to the first agent
              tmp1[fn1.tnfidx()] = tmp2[fn2.tnfidx()];
            if(_PARAM(_IL10dynamics))
              tmp1[fn1.il10idx()] = tmp2[fn2.il10idx()];
          }
          const double errmax = max(abs(err1/yscal1).max(), abs(err2/yscal2).max()) / EPS;
          if(errmax > 1.0)
            {
              double htmp = SAFETY*h*pow(errmax, PSHRINK);
              h = (h >= 0.0 ? max(htmp, 0.1*h) : min(htmp, 0.1*h));
              {
                //tmp = i;
                tmp1 = _initvector;
                tmp2 = other->_initvector;
              }
            }
          else
            {
              if(errmax > ERRCON) dt = SAFETY*h*pow(errmax, PGROW);
              else dt = 5.0*h;
              t += h;
              {
                //i = tmp;
                _initvector = tmp1;
                other->_initvector = tmp2;
              }
              break;
            }
        }
      if((t - t2)*(t2 - t1) >= 0.0)
        {
          writeMembersFromValarray(grid, _initvector);
          other->writeMembersFromValarray(grid, other->_initvector);
          return;
        }
    }
  //Warning! - required too many steps, increase max_steps, safety, or decrease pgrow
  assert(!"Took too many steps");
}
#endif
