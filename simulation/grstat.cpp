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
	: _nMac(0)
	, _nMacResting(0)
	, _nMacInfected(0)
	, _nMacCInfected(0)
	, _nMacActive(0)
	, _nMacDead(0)
	, _nTgam(0)
	, _nTgamActive(0)
	, _nTgamDownRegulated(0)
	, _nTgamDead(0)
	, _nTcyt(0)
	, _nTcytActive(0)
	, _nTcytDownRegulated(0)
	, _nTcytDead(0)
	, _nTreg(0)
	, _nTregActive(0)
	, _nTregDead(0)
	, _totExtMtb(0)
	, _totNonRepExtMtb(0)
	, _totIntMtb(0)
	, _totMacAttractant(0)
	, _totTNF(0)
    , _totIL10(0)
	, _totCCL2(0)
	, _totCCL5(0)
	, _totCXCL9(0)
	, _nApoptosisFasFasL(0)
	, _nMacApoptosisTNF()
	, _nTcellApoptosisTNF(0)
	, _nRestingMacActivationTNF(0)
	, _nInfMacActivationTNF(0)
	, _nMacNFkB(0)
	, _nMacNFkBResting(0)
	, _nMacNFkBInfected(0)
	, _nMacNFkBCInfected(0)
	, _nMacNFkBActive(0)
	, _nMacNFkBDead(0)
	, _nSourceMac(0)
	, _nSourceTgam(0)
	, _nSourceTcyt(0)
	, _nSourceTreg(0)
	, _nMacStat1(0)
	, _nMacStat1Resting(0)
	, _nMacStat1Infected(0)
	, _nMacStat1CInfected(0)
	, _nMacStat1Active(0)
	, _nMacStat1Dead(0)
	, _nMacDeact(0)
	, _nMacDeactResting(0)
	, _nMacDeactInfected(0)
	, _nMacDeactCInfected(0)
	, _nMacDeactActive(0)
	, _nMacDeactDead(0)
	, _nBactAct(0)
	, _areaTNF(0)
	, _areaCellDensity(0)
	, _nCaseated(0)
	, _queueTgam(0)
	, _queueTcyt(0)
	, _queueTreg(0)
	, _fluxTgam(0)
	, _fluxTcyt(0)
	, _fluxTreg(0)
	, _nSourceMacActive(0)
	, _nSourceTgamActive(0)
	, _nSourceTcytActive(0)
	, _nSourceTregActive(0)
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
  for(int i=0;i<NMAC_STATES;i++)
    _macIntMtbStats[i] = o._macIntMtbStats[i];
}
GrStat& GrStat::operator=(const GrStat& o)
{
  memcpy(this, &o, sizeof(GrStat)); //Copy the field values
  _intMtbFreq = new unsigned[_intMtbFreqSize];
  memcpy(_intMtbFreq, o._intMtbFreq, sizeof(unsigned)*_intMtbFreqSize);
  for(int i=0;i<NMAC_STATES;i++)
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
  switch(a->getAgentType())
  {
    case MAC:
    {
      Mac* pMac = static_cast<Mac*>(a);
      assert(pMac != NULL);
      _macIntMtbStats[pMac->getState()](pMac->getIntMtb());
      if(pMac->getState() == MAC_INFECTED || pMac->getState() == MAC_CINFECTED){
        assert(pMac->getIntMtb() < _PARAM(PARAM_MAC_THRESHOLD_BURST_CI_INTMTB));
        assert(pMac->getIntMtb() < _intMtbFreqSize);

        _intMtbFreq[int(pMac->getIntMtb())]++;
      }
      updateMacStatistics(static_cast<Mac*>(pMac)->getState());
      break;
    }
    case TGAM:
      updateTgamStatistics(static_cast<Tgam*>(a)->getState());
      break;
    case TCYT:
      updateTcytStatistics(static_cast<Tcyt*>(a)->getState());
      break;
    case TREG:
      updateTregStatistics(static_cast<Treg*>(a)->getState());
      break;
    default:
      throw std::runtime_error("Unknown agent type");
  }
}
void GrStat::updateMacStatistics(MacState state)
{
	switch (state)
	{
	case MAC_RESTING:
		_nMacResting++;
		break;
	case MAC_INFECTED:
		_nMacInfected++;
		break;
	case MAC_CINFECTED:
		_nMacCInfected++;
		break;
	case MAC_ACTIVE:
		_nMacActive++;
		break;
	case MAC_DEAD:
		_nMacDead++;
		break;
  default: throw std::runtime_error("Unknown Mac state"); break;
	}

	_nMac++;
}

void GrStat::updateMacNFkBStatistics(MacState state)
{
	switch (state)
	{
	case MAC_RESTING:
		_nMacNFkBResting++;
		break;
	case MAC_INFECTED:
		_nMacNFkBInfected++;
		break;
	case MAC_CINFECTED:
		_nMacNFkBCInfected++;
		break;
	case MAC_ACTIVE:
		_nMacNFkBActive++;
		break;
	case MAC_DEAD:
		_nMacNFkBDead++;
		break;
  default: throw std::runtime_error("Unknown Mac state"); break;
	}

	_nMacNFkB++;
}

void GrStat::updateMacStat1Statistics(MacState state)
{
	switch (state)
	{
	case MAC_RESTING:
		_nMacStat1Resting++;
		break;
	case MAC_INFECTED:
		_nMacStat1Infected++;
		break;
	case MAC_CINFECTED:
		_nMacStat1CInfected++;
		break;
	case MAC_ACTIVE:
		_nMacStat1Active++;
		break;
	case MAC_DEAD:
		_nMacStat1Dead++;
		break;
  default: throw std::runtime_error("Unknown Mac state"); break;
	}

	_nMacStat1++;
}

void GrStat::updateMacDeactStatistics(MacState state)
{
	switch (state)
	{
	case MAC_RESTING:
		_nMacDeactResting++;
		break;
	case MAC_INFECTED:
		_nMacDeactInfected++;
		break;
	case MAC_CINFECTED:
		_nMacDeactCInfected++;
		break;
	case MAC_ACTIVE:
		_nMacDeactActive++;
		break;
	case MAC_DEAD:
		_nMacDeactDead++;
		break;
  default: throw std::runtime_error("Unknown Mac state"); break;
	}

	_nMacDeact++;
}

void GrStat::updateTgamStatistics(TgamState state)
{
	switch (state)
	{
	case TGAM_ACTIVE:
		_nTgamActive++;
		break;
	case TGAM_DOWN_REGULATED:
		_nTgamDownRegulated++;
		break;
	case TGAM_DEAD:
		_nTgamDead++;
		break;
  default: throw std::runtime_error("Unknown Tgam state"); break;
	}

	_nTgam++;
}

void GrStat::updateTcytStatistics(TcytState state)
{
	switch (state)
	{
	case TCYT_ACTIVE:
		_nTcytActive++;
		break;
	case TCYT_DOWN_REGULATED:
		_nTcytDownRegulated++;
		break;
	case TCYT_DEAD:
		_nTcytDead++;
		break;
  default: throw std::runtime_error("Unknown Tcyt state"); break;
	}

	_nTcyt++;
}

void GrStat::updateTregStatistics(TregState state)
{
	switch (state)
	{
	case TREG_ACTIVE:
		_nTregActive++;
		break;
	case TREG_DEAD:
		_nTregDead++;
		break;
  default: throw std::runtime_error("Unknown Treg state"); break;
	}

	_nTreg++;
}

void GrStat::resetAgentStats()
{
	_nMac = _nMacResting = _nMacInfected = 
		_nMacCInfected = _nMacDead = _nMacActive = 0;
	
	_nMacNFkB = _nMacNFkBResting = _nMacNFkBInfected =
		_nMacNFkBCInfected = _nMacNFkBDead = _nMacNFkBActive = 0;
	
	_nMacStat1 = _nMacStat1Resting = _nMacStat1Infected = _nMacStat1CInfected = 
		_nMacStat1Dead = _nMacStat1Active = 0;
	
	_nMacDeact = _nMacDeactResting = _nMacDeactInfected = _nMacDeactCInfected =
		_nMacDeactDead = _nMacDeactActive = 0;

	_nTgam = _nTgamActive = _nTgamDead = _nTgamDownRegulated = 0;
	
	_nTcyt = _nTcytDead = _nTcytDownRegulated = _nTcytActive = 0;
	
	_nTreg = _nTregDead = _nTregActive = 0;

  memset(_nMacApoptosisTNF, 0 , sizeof(int)*NMAC_STATES);
  memset(_intMtbFreq, 0, sizeof(unsigned)*_intMtbFreqSize);
  for(int i=0;i<NMAC_STATES;i++)
    _macIntMtbStats[i] = Stat();
}

void GrStat::reset()
{
	_totIntMtb = _totExtMtb = _totNonRepExtMtb = _totMacAttractant = _totTNF = _totCCL2 = _totCCL5 = _totCXCL9 = _totIL10 = 0;
	
	_nSourceMac = _nSourceTcyt = _nSourceTgam = _nSourceTreg = 0;
	
	_nSourceMacActive = _nSourceTcytActive = _nSourceTgamActive = _nSourceTregActive = 0;

	_nCaseated = 0;

	_areaTNF = _areaCellDensity = 0;
}

void GrStat::serialize(std::ostream& out) const
{
	assert(out.good());

	Serialization::writeHeader(out, GrStat::_ClassName);

	out << _nMac << std::endl;
	out << _nMacResting << std::endl;
	out << _nMacInfected << std::endl;
	out << _nMacCInfected << std::endl;
	out << _nMacActive << std::endl;
	out << _nMacDead << std::endl;

	out << _nTgam << std::endl;
	out << _nTgamActive << std::endl;
	out << _nTgamDownRegulated << std::endl;
	out << _nTgamDead << std::endl;

	out << _nTcyt << std::endl;
	out << _nTcytActive << std::endl;
	out << _nTcytDownRegulated << std::endl;
	out << _nTcytDead << std::endl;

	out << _nTreg << std::endl;
	out << _nTregActive << std::endl;
	out << _nTregDead << std::endl;

	out << _totExtMtb << std::endl;
	out << _totNonRepExtMtb << std::endl;
	out << _totIntMtb << std::endl;

	out << _totMacAttractant << std::endl;
	out << _totTNF << std::endl;
    out << _totIL10 << std::endl;
	out << _totCCL2 << std::endl;
	out << _totCCL5 << std::endl;
	out << _totCXCL9 << std::endl;

	out << _nApoptosisFasFasL << std::endl;
  for(unsigned i=0;i<NMAC_STATES;i++)
  	out << _nMacApoptosisTNF[i] << std::endl;

	out << _nTcellApoptosisTNF << std::endl;
	
	out << _nRestingMacActivationTNF << std::endl;
	out << _nInfMacActivationTNF << std::endl;

	out << _nMacNFkB << std::endl;
	out << _nMacNFkBResting << std::endl;
	out << _nMacNFkBInfected << std::endl;
	out << _nMacNFkBCInfected << std::endl;
	out << _nMacNFkBActive << std::endl;
	out << _nMacNFkBDead << std::endl;

	out << _nSourceMac << std::endl;
	out << _nSourceTgam << std::endl;
	out << _nSourceTcyt << std::endl;
	out << _nSourceTreg << std::endl;

	out << _nMacStat1 << std::endl;
	out << _nMacStat1Resting << std::endl;
	out << _nMacStat1Infected << std::endl;
	out << _nMacStat1CInfected << std::endl;
	out << _nMacStat1Active << std::endl;
	out << _nMacStat1Dead << std::endl;

	out << _nMacDeact << std::endl;
	out << _nMacDeactResting << std::endl;
	out << _nMacDeactInfected << std::endl;
	out << _nMacDeactCInfected << std::endl;
	out << _nMacDeactActive << std::endl;
	out << _nMacDeactDead << std::endl;

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
	out << _fluxTgam << std::endl;
	out << _fluxTcyt << std::endl;
	out << _fluxTreg << std::endl;
	out << _nSourceMacActive << std::endl;
	out << _nSourceTgamActive << std::endl;
	out << _nSourceTcytActive << std::endl;
	out << _nSourceTregActive << std::endl;
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

	Serialization::writeFooter(out, GrStat::_ClassName);
}

void GrStat::deserialize(std::istream& in)
{
	assert(in.good());

	if (!Serialization::readHeader(in, GrStat::_ClassName))
	{
		exit(1);
	}

	in >>_nMac;
	in >>_nMacResting;
	in >>_nMacInfected;
	in >>_nMacCInfected;
	in >>_nMacActive;
	in >>_nMacDead;

	in >>_nTgam;
	in >>_nTgamActive;
	in >>_nTgamDownRegulated;
	in >>_nTgamDead;

	in >>_nTcyt;
	in >>_nTcytActive;
	in >>_nTcytDownRegulated;
	in >>_nTcytDead;

	in >>_nTreg;
	in >>_nTregActive;
	in >>_nTregDead;

	in >>_totExtMtb;
	in >>_totNonRepExtMtb;
	in >>_totIntMtb;

	in >>_totMacAttractant;
	in >>_totTNF;
    in >>_totIL10;
	in >>_totCCL2;
	in >>_totCCL5;
	in >>_totCXCL9;

	in >>_nApoptosisFasFasL;
  for(unsigned i=0;i<NMAC_STATES;i++)
  	in >>_nMacApoptosisTNF[i];
	
  in >> _nTcellApoptosisTNF;
	
	in >> _nRestingMacActivationTNF;
	in >> _nInfMacActivationTNF;

	in >>_nMacNFkB;
	in >>_nMacNFkBResting;
	in >>_nMacNFkBInfected;
	in >>_nMacNFkBCInfected;
	in >>_nMacNFkBActive;
	in >>_nMacNFkBDead;

	in >>_nSourceMac;
	in >>_nSourceTgam;
	in >>_nSourceTcyt;
	in >>_nSourceTreg;

	in >>_nMacStat1;
	in >>_nMacStat1Resting;
	in >>_nMacStat1Infected;
	in >>_nMacStat1CInfected;
	in >>_nMacStat1Active;
	in >>_nMacStat1Dead;

	in >>_nMacDeact;
	in >>_nMacDeactResting;
	in >>_nMacDeactInfected;
	in >>_nMacDeactCInfected;
	in >>_nMacDeactActive;
	in >>_nMacDeactDead;

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
	in >>_fluxTgam;
	in >>_fluxTcyt;
	in >>_fluxTreg;
	in >>_nSourceMacActive;
	in >>_nSourceTgamActive;
	in >>_nSourceTcytActive;
	in >>_nSourceTregActive;
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

	if (!Serialization::readFooter(in, GrStat::_ClassName))
	{
		exit(1);
	}
}
