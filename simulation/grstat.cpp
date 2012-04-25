/*
 * grstat.cpp
 *
 *  Created on: 30-nov-2009
 *      Author: M. El-Kebir
 */

#include "grstat.h"
#include "params.h"
#include "serialization.h"
#include "macrophage.h"
#include "tgamma.h"
#include "tcytotoxic.h"
#include "tregulatory.h"

const std::string GrStat::_ClassName = "GrStat";

GrStat::GrStat()
	: _nAgents()
	, _nMac()
	, _nTgam()
	, _nTcyt()
	, _nTreg()
	, _totExtMtb(0)
	, _totNonRepExtMtb(0)
	, _totIntMtb(0)
	, _totMacAttractant(0)
	, _totTNF(0)
    , _totIL10(0)
	, _totCCL2(0)
	, _totCCL5(0)
	, _totCXCL9(0)
    , _totTNFR1int(0)
    , _totkmRNA(0)
	, _nApoptosisFasFasL(0)
	, _nMacApoptosisTNF()
	, _nTcellApoptosisTNF(0)
	, _nRestingMacActivationTNF(0)
	, _nInfMacActivationTNF(0)
	, _nMacNFkB()
	, _nSource()
	, _nSourceActive()
	, _nSourceMacCrowded(0)
	, _nSourceTgamCrowded(0)
	, _nSourceTcytCrowded(0)
	, _nSourceTregCrowded(0)
	, _nMacStat1()
	, _nMacDeact()
	, _nBactAct(0)
	, _areaTNF(0)
	, _areaCellDensity(0)
	, _nCaseated(0)
	, _grStatus()
	, _queueTgam(0)
	, _queueTcyt(0)
	, _queueTreg(0)
	, _queueTgamDie(0)
	, _queueTcytDie(0)
	, _queueTregDie(0)
	, _recruitedTgam(0)
	, _recruitedTcyt(0)
	, _recruitedTreg(0)
	, _fluxTgam(0)
	, _fluxTcyt(0)
	, _fluxTreg(0)
	, _MDC(0)
	, _N4(0)
	, _TH0(0)
	, _TH1(0)
	, _N8(0)
	, _T80(0)
	, _T8(0)
	, _TC(0)
	, _TH0lung(0)
	, _TH1lung(0)
	, _T80lung(0)
	, _T8lung(0)
	, _TClung(0)
	, _nCellTnfInhibit(0)

		// + 1 because burst threshold might not be an integer value, and if not,
		// when a macrophage's intMtb count is truncated to the floor(burst threshold)
		// and then used as an index into _intMtbFreq an array over reference will occur.
		// Ex. PARAM_MAC_THRESHOLD_BURST_CI_INTMTB 22.4
		// If a mac's intMtb grows to 22.1, when it's intMtb value, as the integer 22, is used as an index
		// into _intMtbFreq, so _intMtbFreq must be dimensioned to 23 or an array over reference will occur.
	, _intMtbFreqSize(int(_PARAM(PARAM_MAC_THRESHOLD_BURST_CI_INTMTB)) + 1)
	, _intMtbFreq(new unsigned[_intMtbFreqSize])
{
	for (int i = 0; i < NOUTCOMES; i++)
	{
		_grStatus[i] = GR_NONE;
	}

	memset(_intMtbFreq, 0, sizeof(unsigned)*_intMtbFreqSize);
}
GrStat::GrStat(const GrStat& o)
{
  memcpy(this, &o, sizeof(GrStat)); //Copy the field values
  _intMtbFreq = new unsigned[_intMtbFreqSize];
  memcpy(_intMtbFreq, o._intMtbFreq, sizeof(unsigned)*_intMtbFreqSize);
  for(int i=0;i<Mac::NSTATES;i++)
    _macIntMtbStats[i] = o._macIntMtbStats[i];
}
GrStat& GrStat::operator=(const GrStat& o)
{
  memcpy(this, &o, sizeof(GrStat)); //Copy the field values
  _intMtbFreq = new unsigned[_intMtbFreqSize];
  memcpy(_intMtbFreq, o._intMtbFreq, sizeof(unsigned)*_intMtbFreqSize);
  for(int i=0;i<Mac::NSTATES;i++)
    _macIntMtbStats[i] = o._macIntMtbStats[i];
  return *this;
}

GrStat::~GrStat()
{
  if(_intMtbFreq)
    delete[] _intMtbFreq;
}

void GrStat::updateAgentStatistics(Agent* a)
{
  assert(a!=NULL);
  ++_nAgents[a->getAgentType()];
  switch(a->getAgentType())
  {
    case MAC:
    {
      Mac* pMac = static_cast<Mac*>(a);
      assert(pMac != NULL);
      _macIntMtbStats[pMac->getState()](pMac->getIntMtb());
      if(pMac->getState() == Mac::MAC_INFECTED || pMac->getState() == Mac::MAC_CINFECTED){

    	// This can happen if  _PARAM(PARAM_MAC_THRESHOLD_BURST_CI_INTMTB) < _PARAM(PARAM_MAC_THRESHOLD_BECOME_CI_INTMTB)
    	// or if the intMtb growth rate is high enough for intMtb for a mac > both PARAM(PARAM_MAC_THRESHOLD_BECOME_CI_INTMTB)
    	// and _PARAM(PARAM_MAC_THRESHOLD_BURST_CI_INTMTB). In either case intMtb >  _PARAM(PARAM_MAC_THRESHOLD_BURST_CI_INTMTB)
    	// when an Mi becomes an MCi. The new MCi won't burst until the next time step, so the assert condition would be true here.
        //assert(pMac->getIntMtb() < _PARAM(PARAM_MAC_THRESHOLD_BURST_CI_INTMTB));

        assert(pMac->getIntMtb() < _intMtbFreqSize);

        _intMtbFreq[int(pMac->getIntMtb())]++;
      }
      ++_nMac[a->getState()];
      break;
    }
    case TGAM:
      ++_nTgam[a->getState()];
      break;
    case TCYT:
      ++_nTcyt[a->getState()];
      break;
    case TREG:
      ++_nTreg[a->getState()];
      break;
    default:
      throw std::runtime_error("Unknown agent type");
  }
}

void GrStat::updateMacNFkBStatistics(Mac::State state)
{
  assert(state >= 0 && state < Mac::NSTATES);
  ++_nMacNFkB[state];
}

void GrStat::updateMacStat1Statistics(Mac::State state)
{
  assert(state >= 0 && state < Mac::NSTATES);
  ++_nMacStat1[state];
}

void GrStat::updateMacDeactStatistics(Mac::State state)
{
  assert(state >= 0 && state < Mac::NSTATES);
  ++_nMacDeact[state];
}

void GrStat::resetAgentStats()
{
  memset(_nAgents, 0, sizeof(int)*NAGENTS);
  memset(_nMac, 0, sizeof(int)*Mac::NSTATES);
  memset(_nMacNFkB, 0, sizeof(int)*Mac::NSTATES);
  memset(_nMacStat1, 0, sizeof(int)*Mac::NSTATES);
  memset(_nMacDeact, 0, sizeof(int)*Mac::NSTATES);
  memset(_nTgam, 0, sizeof(int)*Tgam::NSTATES);
  memset(_nTcyt, 0, sizeof(int)*Tcyt::NSTATES);
  memset(_nTreg, 0, sizeof(int)*Treg::NSTATES);

  memset(_nMacApoptosisTNF, 0 , sizeof(int)*Mac::NSTATES);
  memset(_intMtbFreq, 0, sizeof(unsigned)*_intMtbFreqSize);
  for(int i=0;i<Mac::NSTATES;i++)
    _macIntMtbStats[i] = Stat();
}

void GrStat::reset()
{
	_totIntMtb = _totExtMtb = _totNonRepExtMtb = _totMacAttractant = _totTNF = _totCCL2 = _totCCL5 = _totCXCL9 = _totIL10 = 0;
	
  memset(_nSource, 0, sizeof(int)*NAGENTS);
  memset(_nSourceActive, 0, sizeof(int)*NAGENTS);
	
	_nCaseated = 0;

	_queueTgam = _queueTcyt = _queueTreg = 0;

	_queueTgamDie = _queueTcytDie = _queueTregDie = 0;

	_recruitedTgam = _recruitedTcyt = _recruitedTreg = 0;

	_areaTNF = _areaCellDensity = 0;
	
	_nCellTnfInhibit = _totTNFR1int = _totkmRNA = 0;
}

void GrStat::serialize(std::ostream& out) const
{
	assert(out.good());

	Serialization::writeHeader(out, GrStat::_ClassName);

	out <<_nAgents[MAC] << std::endl;
  for(size_t i=0;i<Mac::NSTATES;i++)
    out <<_nMac[i] << std::endl;

	out <<_nAgents[TGAM] << std::endl;
  for(size_t i=0;i<Tgam::NSTATES;i++)
    out <<_nTgam[i] << std::endl;

	out <<_nAgents[TCYT] << std::endl;
  for(size_t i=0;i<Tcyt::NSTATES;i++)
    out <<_nTcyt[i] << std::endl;

	out <<_nAgents[TREG] << std::endl;
  for(size_t i=0;i<Treg::NSTATES;i++)
    out <<_nTreg[i] << std::endl;

	out << _totExtMtb << std::endl;
	out << _totNonRepExtMtb << std::endl;
	out << _totIntMtb << std::endl;

	out << _totMacAttractant << std::endl;
	out << _totTNF << std::endl;
    out << _totIL10 << std::endl;
	out << _totCCL2 << std::endl;
	out << _totCCL5 << std::endl;
	out << _totCXCL9 << std::endl;
    out << _totTNFR1int << std::endl;
    out << _totkmRNA << std::endl;

	out << _nApoptosisFasFasL << std::endl;
	
  for(unsigned i=0;i<Mac::NSTATES;i++)
  	out << _nMacApoptosisTNF[i] << std::endl;

	out << _nTcellApoptosisTNF << std::endl;
	
	out << _nRestingMacActivationTNF << std::endl;
	out << _nInfMacActivationTNF << std::endl;

  for(size_t i=0;i<Mac::NSTATES;i++)
    out <<_nMacNFkB[i] << std::endl;

  for(size_t i=0;i<NAGENTS;i++)
    out <<_nSource[i] << std::endl;

  for(size_t i=0;i<NAGENTS;i++)
    out << _nSourceActive[i] << std::endl;

  out << _nSourceMacCrowded << std::endl;
  out << _nSourceTgamCrowded << std::endl;
  out << _nSourceTcytCrowded << std::endl;
  out << _nSourceTregCrowded << std::endl;

  for(size_t i=0;i<Mac::NSTATES;i++)
    out <<_nMacStat1[i] << std::endl;

  for(size_t i=0;i<Mac::NSTATES;i++)
    out <<_nMacDeact[i] << std::endl;

	out << _nBactAct << std::endl;
	out << _areaTNF << std::endl;
	out << _areaCellDensity << std::endl;
	out << _nCaseated << std::endl;

	for (int i = 0; i < NOUTCOMES; i++)
	{
		int intVal = (int) _grStatus[i];
		out << intVal << std::endl;
	}

	out << _queueTgam << std::endl;
	out << _queueTcyt << std::endl;
	out << _queueTreg << std::endl;

	out << _queueTgamDie << std::endl;
	out << _queueTcytDie << std::endl;
	out << _queueTregDie << std::endl;

	out << _recruitedTgam << std::endl;
	out << _recruitedTcyt << std::endl;
	out << _recruitedTreg << std::endl;

	out << _fluxTgam << std::endl;
	out << _fluxTcyt << std::endl;
	out << _fluxTreg << std::endl;

	out << _MDC << std::endl;
	out << _N4 << std::endl;
	out << _TH0 << std::endl;
	out << _TH1 << std::endl;
	out << _N8 << std::endl;
	out << _T80 << std::endl;
	out << _T8 << std::endl;
	out << _TC << std::endl;
	out << _TH0lung << std::endl;
	out << _TH1lung << std::endl;
	out << _T80lung << std::endl;
	out << _T8lung << std::endl;
	out << _TClung << std::endl;
	out << _nCellTnfInhibit << std::endl;

	out << _intMtbFreqSize << std::endl;

	for (size_t i = 0; i < _intMtbFreqSize; i++)
	{
		out << _intMtbFreq[i] << std::endl;
	}

	Serialization::writeFooter(out, GrStat::_ClassName);
}

void GrStat::deserialize(std::istream& in)
{
	assert(in.good());

	if (!Serialization::readHeader(in, GrStat::_ClassName))
	{
		exit(1);
	}

	in >>_nAgents[MAC];
  for(size_t i=0;i<Mac::NSTATES;i++)
    in >>_nMac[i];

	in >>_nAgents[TGAM];
  for(size_t i=0;i<Tgam::NSTATES;i++)
    in >>_nTgam[i];

	in >>_nAgents[TCYT];
  for(size_t i=0;i<Tcyt::NSTATES;i++)
    in >>_nTcyt[i];

	in >>_nAgents[TREG];
  for(size_t i=0;i<Treg::NSTATES;i++)
    in >>_nTreg[i];

	in >>_totExtMtb;
	in >>_totNonRepExtMtb;
	in >>_totIntMtb;

	in >>_totMacAttractant;
	in >>_totTNF;
    in >>_totIL10;
	in >>_totCCL2;
	in >>_totCCL5;
	in >>_totCXCL9;
    in >>_totTNFR1int;
    in >>_totkmRNA;

	in >>_nApoptosisFasFasL;

  for(unsigned i=0;i<Mac::NSTATES;i++)
  	in >>_nMacApoptosisTNF[i];
	
    in >> _nTcellApoptosisTNF;
	
	in >> _nRestingMacActivationTNF;
	in >> _nInfMacActivationTNF;

  for(size_t i=0;i<Mac::NSTATES;i++)
    in >>_nMacNFkB[i];

  for(size_t i=0;i<NAGENTS;i++)
    in >>_nSource[i];

  for(size_t i=0;i<NAGENTS;i++)
    in >>_nSourceActive[i];

  in >> _nSourceMacCrowded;
  in >> _nSourceTgamCrowded;
  in >> _nSourceTcytCrowded;
  in >> _nSourceTregCrowded;

  for(size_t i=0;i<Mac::NSTATES;i++)
    in >>_nMacStat1[i];

  for(size_t i=0;i<Mac::NSTATES;i++)
    in >>_nMacDeact[i];

	in >>_nBactAct;
	in >>_areaTNF;
	in >>_areaCellDensity;
	in >>_nCaseated;


	for (int i = 0; i < NOUTCOMES; i++)
	{
		int intVal;
		in >> intVal;
		_grStatus[i] = (GrStatus) intVal;
	}

	in >>_queueTgam;
	in >>_queueTcyt;
	in >>_queueTreg;

	in >> _queueTgamDie;
	in >> _queueTcytDie;
	in >> _queueTregDie;

	in >> _recruitedTgam;
	in >> _recruitedTcyt;
	in >> _recruitedTreg;

	in >>_fluxTgam;
	in >>_fluxTcyt;
	in >>_fluxTreg;

	in >>_MDC;
	in >>_N4;
	in >>_TH0;
	in >>_TH1;
	in >>_N8;
	in >>_T80;
	in >>_T8;
	in >>_TC;
	in >>_TH0lung;
	in >>_TH1lung;
	in >>_T80lung;
	in >>_T8lung;
	in >>_TClung;
	in >> _nCellTnfInhibit;

	in >> _intMtbFreqSize;

	for (size_t i = 0; i < _intMtbFreqSize; i++)
	{
		in >> _intMtbFreq[i];
	}

	if (!Serialization::readFooter(in, GrStat::_ClassName))
	{
		exit(1);
	}
}
