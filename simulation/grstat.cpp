/*
 * grstat.cpp
 *
 *  Created on: 30-nov-2009
 *      Author: M. El-Kebir
 */

#include "grstat.h"

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
	, _totIntMtb(0)
	, _totMacAttractant(0)
	, _totTNF(0)
	, _totCCL2(0)
	, _totCCL5(0)
	, _totCXCL9(0)
	, _nApoptosisFasFasL(0)
	, _nApoptosisTNF(0)
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
	, _area(0)
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
{
	for (int i = 0; i < NOUTCOMES; i++)
	{
		_grStatus[i] = GR_NONE;
	}
}

GrStat::~GrStat()
{
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
}

void GrStat::reset()
{
	_totIntMtb = _totExtMtb = _totMacAttractant = _totTNF = _totCCL2 = _totCCL5 = _totCXCL9 = 0;
	
	_nSourceMac = _nSourceTcyt = _nSourceTgam = _nSourceTreg = 0;
	
	_nSourceMacActive = _nSourceTcytActive = _nSourceTgamActive = _nSourceTregActive = 0;

	_area = 0;
}

void GrStat::serialize(std::ostream& out) const
{
	assert(out.good());

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
	out << _totIntMtb << std::endl;
	out << _totMacAttractant << std::endl;
	out << _totTNF << std::endl;
	out << _totCCL2 << std::endl;
	out << _totCCL5 << std::endl;
	out << _totCXCL9 << std::endl;
	out << _nApoptosisFasFasL << std::endl;
	out << _nApoptosisTNF << std::endl;
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
	out << _nBactAct << std::endl;
	out << _area << std::endl;

	for (int i = 0; i < NOUTCOMES; i++)
	{
		int intVal = (int) _grStatus[i];
		out << intVal << std::endl;
	}
}

void GrStat::deserialize(std::istream& in)
{
	assert(in.good());

	in >> _nMac;
	in >> _nMacResting;
	in >> _nMacInfected;
	in >> _nMacCInfected;
	in >> _nMacActive;
	in >> _nMacDead;
	in >> _nTgam;
	in >> _nTgamActive;
	in >> _nTgamDownRegulated;
	in >> _nTgamDead;
	in >> _nTcyt;
	in >> _nTcytActive;
	in >> _nTcytDownRegulated;
	in >> _nTcytDead;
	in >> _nTreg;
	in >> _nTregActive;
	in >> _nTregDead;
	in >> _totExtMtb;
	in >> _totIntMtb;
	in >> _totMacAttractant;
	in >> _totTNF;
	in >> _totCCL2;
	in >> _totCCL5;
	in >> _totCXCL9;
	in >> _nApoptosisFasFasL;
	in >> _nApoptosisTNF;
	in >> _nMacNFkB;
	in >> _nMacNFkBResting;
	in >> _nMacNFkBInfected;
	in >> _nMacNFkBCInfected;
	in >> _nMacNFkBActive;
	in >> _nMacNFkBDead;
	in >> _nSourceMac;
	in >> _nSourceTgam;
	in >> _nSourceTcyt;
	in >> _nSourceTreg;
	in >> _nMacStat1;
	in >> _nMacStat1Resting;
	in >> _nMacStat1Infected;
	in >> _nMacStat1CInfected;
	in >> _nMacStat1Active;
	in >> _nMacStat1Dead;
	in >> _nBactAct;
	in >> _area;

	for (int i = 0; i < NOUTCOMES; i++)
	{
		int intVal;
		in >> intVal;
		_grStatus[i] = (GrStatus) intVal;
	}
}
