/*
 * grstat.h
 *
 *  Created on: 30-nov-2009
 *      Author: M. El-Kebir
 */

#ifndef GRSTAT_H
#define GRSTAT_H

#include "gr.h"

typedef enum {	GR_CONTAINMENT, GR_CONTAINMENT_INCONSISTENT, GR_CLEARANCE,
				GR_DISSEMINATION, GR_DISSEMINATION_INCONSISTENT, GR_UNKNOWN, GR_NONE} GrStatus;

class GrStat
{
private:
	int _nMac;
	int _nMacResting;
	int _nMacInfected;
	int _nMacCInfected;
	int _nMacActive;
	int _nMacDead;
	int _nTgam;
	int _nTgamActive;
	int _nTgamDownRegulated;
	int _nTgamDead;
	int _nTcyt;
	int _nTcytActive;
	int _nTcytDownRegulated;
	int _nTcytDead;
	int _nTreg;
	int _nTregActive;
	int _nTregDead;
	double _totExtMtb;
	double _totIntMtb;
	double _totMacAttractant;
	double _totTNF;
	double _totCCL2;
	double _totCCL5;
	double _totCXCL9;
	int _nApoptosisFasFasL;
	int _nApoptosisTNF;
	int _nMacNFkB;
	int _nMacNFkBResting;
	int _nMacNFkBInfected;
	int _nMacNFkBCInfected;
	int _nMacNFkBActive;
	int _nMacNFkBDead;
	int _nSourceMac;
	int _nSourceTgam;
	int _nSourceTcyt;
	int _nSourceTreg;
	int _nMacStat1;
	int _nMacStat1Resting;
	int _nMacStat1Infected;
	int _nMacStat1CInfected;
	int _nMacStat1Active;
	int _nMacStat1Dead;
	int _nMacDeact;
	int _nMacDeactResting;
	int _nMacDeactInfected;
	int _nMacDeactCInfected;
	int _nMacDeactActive;
	int _nMacDeactDead;
	int _nBactAct;
	int _area;
	GrStatus _grStatus[NOUTCOMES];

public:
	GrStat();
	~GrStat();
	int getNrOfMac() const;
	int getNrOfTgam() const;
	int getNrOfTcyt() const;
	int getNrOfTreg() const;
	int getNrOfMacResting() const;
	int getNrOfMacInfected() const;
	int getNrOfMacCInfected() const;
	int getNrOfMacActive() const;
	int getNrOfMacDead() const;
	int getNrOfMacNFkB() const;
	int getNrOfMacNFkBResting() const;
	int getNrOfMacNFkBInfected() const;
	int getNrOfMacNFkBCInfected() const;
	int getNrOfMacNFkBActive() const;
	int getNrOfMacNFkBDead() const;
	int getNrOfMacStat1() const;
	int getNrOfMacStat1Resting() const;
	int getNrOfMacStat1Infected() const;
	int getNrOfMacStat1CInfected() const;
	int getNrOfMacStat1Active() const;
	int getNrOfMacStat1Dead() const;
	int getNrOfMacDeact() const;
	int getNrOfMacDeactResting() const;
	int getNrOfMacDeactInfected() const;
	int getNrOfMacDeactCInfected() const;
	int getNrOfMacDeactActive() const;
	int getNrOfMacDeactDead() const;
	int getNrOfTgamActive() const;
	int getNrOfTgamDownRegulated() const;
	int getNrOfTgamDead() const;
	int getNrOfTcytActive() const;
	int getNrOfTcytDead() const;
	int getNrOfTcytDownRegulated() const;
	int getNrOfTregActive() const;
	int getNrOfTregDead() const;
	double getTotExtMtb() const;
	double getTotIntMtb() const;
	double getTotMacAttractant() const;
	double getTotTNF() const;
	double getTotCCL2() const;
	double getTotCCL5() const;
	double getTotCXCL9() const;
	void updateMacStatistics(MacState state);
	void updateTgamStatistics(TgamState state);
	void updateTcytStatistics(TcytState state);
	void updateTregStatistics(TregState state);
	void updateMacNFkBStatistics(MacState state);
	void updateMacStat1Statistics(MacState state);
	void updateMacDeactStatistics(MacState state);
	void reset();
	void incTotExtMtb(double dExtMtb);
	void incTotIntMtb(double dIntMtb);
	void incTotMacAttractant(double dMacAttractant);
	void incTotTNF(double dTNF);
	void incTotCCL2(double dCCL2);
	void incTotCCL5(double dCCL5);
	void incTotCXCL9(double dCXCL9);
	int getNrApoptosisFasFasL() const;
	int getNrApoptosisTNF() const;
	void incApoptosisFasFasL();
	void incApoptosisTNF();
	int getNrSourcesMac() const;
	int getNrSourcesTgam() const;
	int getNrSourcesTcyt() const;
	int getNrSourcesTreg() const;
	void incNrSourcesMac();
	void incNrSourcesTgam();
	void incNrSourcesTcyt();
	void incNrSourcesTreg();
	int getNrBactAct() const;
	void incNrBactAct();
	int getArea() const;
	void incArea();
	GrStatus getGrStatus(int index) const;
	void setGrStatus(int index, GrStatus status);
	void serialize(std::ostream& out) const;
	void deserialize(std::istream& in);
};

inline int GrStat::getNrBactAct() const
{
	return _nBactAct;
}

inline void GrStat::incNrBactAct()
{
	_nBactAct++;
}

inline GrStatus GrStat::getGrStatus(int index) const
{
	assert(0 <= index && index < NOUTCOMES);
	return _grStatus[index];
}

inline void GrStat::setGrStatus(int index, GrStatus status)
{
	assert(0 <= index && index < NOUTCOMES);
	_grStatus[index] = status;
}

inline int GrStat::getArea() const
{
	return _area;
}

inline void GrStat::incArea()
{
	_area++;
}

inline int GrStat::getNrSourcesMac() const
{
	return _nSourceMac;
}

inline int GrStat::getNrSourcesTgam() const
{
	return _nSourceTgam;
}

inline int GrStat::getNrSourcesTcyt() const
{
	return _nSourceTcyt;
}

inline int GrStat::getNrSourcesTreg() const
{
	return _nSourceTreg;
}

inline void GrStat::incNrSourcesMac()
{
	_nSourceMac++;
}

inline void GrStat::incNrSourcesTgam()
{
	_nSourceTgam++;
}

inline void GrStat::incNrSourcesTcyt()
{
	_nSourceTcyt++;
}

inline void GrStat::incNrSourcesTreg()
{
	_nSourceTreg++;
}

inline int GrStat::getNrApoptosisFasFasL() const
{
	return _nApoptosisFasFasL;
}

inline int GrStat::getNrApoptosisTNF() const
{
	return _nApoptosisTNF;
}

inline void GrStat::incApoptosisFasFasL()
{
	_nApoptosisFasFasL++;
}

inline void GrStat::incApoptosisTNF()
{
	_nApoptosisTNF++;
}

inline void GrStat::incTotExtMtb(double dExtMtb)
{
	_totExtMtb += dExtMtb;
}

inline void GrStat::incTotIntMtb(double dIntMtb)
{
	_totIntMtb += dIntMtb;
}

inline void GrStat::incTotMacAttractant(double dMacAttractant)
{
	_totMacAttractant += dMacAttractant;
}

inline void GrStat::incTotTNF(double dTNF)
{
	_totTNF += dTNF;
}

inline void GrStat::incTotCCL2(double dCCL2)
{
	_totCCL2 += dCCL2;
}

inline void GrStat::incTotCCL5(double dCCL5)
{
	_totCCL5 += dCCL5;
}

inline void GrStat::incTotCXCL9(double dCXCL9)
{
	_totCXCL9 += dCXCL9;
}

inline double GrStat::getTotMacAttractant() const
{
	return _totMacAttractant;
}

inline double GrStat::getTotTNF() const
{
	return _totTNF;
}

inline double GrStat::getTotCCL2() const
{
	return _totCCL2;
}

inline double GrStat::getTotCCL5() const
{
	return _totCCL5;
}

inline double GrStat::getTotCXCL9() const
{
	return _totCXCL9;
}

inline double GrStat::getTotExtMtb() const
{
	return _totExtMtb;
}

inline double GrStat::getTotIntMtb() const
{
	return _totIntMtb;
}

inline int GrStat::getNrOfMacResting() const
{
	return _nMacResting;
}

inline int GrStat::getNrOfMacInfected() const
{
	return _nMacInfected;
}

inline int GrStat::getNrOfMacCInfected() const
{
	return _nMacCInfected;
}

inline int GrStat::getNrOfMacActive() const
{
	return _nMacActive;
}

inline int GrStat::getNrOfMacDead() const
{
	return _nMacDead;
}

inline int GrStat::getNrOfMac() const
{
	return _nMac;
}

inline int GrStat::getNrOfMacNFkBResting() const
{
	return _nMacNFkBResting;
}

inline int GrStat::getNrOfMacNFkBInfected() const
{
	return _nMacNFkBInfected;
}

inline int GrStat::getNrOfMacNFkBCInfected() const
{
	return _nMacNFkBCInfected;
}

inline int GrStat::getNrOfMacNFkBActive() const
{
	return _nMacNFkBActive;
}

inline int GrStat::getNrOfMacNFkBDead() const
{
	return _nMacNFkBDead;
}

inline int GrStat::getNrOfMacNFkB() const
{
	return _nMacNFkB;
}

inline int GrStat::getNrOfTgamActive() const
{
	return _nTgamActive;
}

inline int GrStat::getNrOfTgamDead() const
{
	return _nTgamDead;
}

inline int GrStat::getNrOfTgamDownRegulated() const
{
	return _nTgamDownRegulated;
}

inline int GrStat::getNrOfTgam() const
{
	return _nTgam;
}

inline int GrStat::getNrOfTcytDead() const
{
	return _nTcytDead;
}

inline int GrStat::getNrOfTcytActive() const
{
	return _nTcytActive;
}

inline int GrStat::getNrOfTcytDownRegulated() const
{
	return _nTcytDownRegulated;
}

inline int GrStat::getNrOfTcyt() const
{
	return _nTcyt;
}

inline int GrStat::getNrOfTregDead() const
{
	return _nTregDead;
}

inline int GrStat::getNrOfTregActive() const
{
	return _nTregActive;
}

inline int GrStat::getNrOfTreg() const
{
	return _nTreg;
}

inline int GrStat::getNrOfMacStat1Resting() const
{
	return _nMacStat1Resting;
}

inline int GrStat::getNrOfMacStat1Infected() const
{
	return _nMacStat1Infected;
}

inline int GrStat::getNrOfMacStat1CInfected() const
{
	return _nMacStat1CInfected;
}

inline int GrStat::getNrOfMacStat1Active() const
{
	return _nMacStat1Active;
}

inline int GrStat::getNrOfMacStat1Dead() const
{
	return _nMacStat1Dead;
}

inline int GrStat::getNrOfMacStat1() const
{
	return _nMacStat1;
}

inline int GrStat::getNrOfMacDeactResting() const
{
	return _nMacDeactResting;
}

inline int GrStat::getNrOfMacDeactInfected() const
{
	return _nMacDeactInfected;
}

inline int GrStat::getNrOfMacDeactCInfected() const
{
	return _nMacDeactCInfected;
}

inline int GrStat::getNrOfMacDeactActive() const
{
	return _nMacDeactActive;
}

inline int GrStat::getNrOfMacDeactDead() const
{
	return _nMacDeactDead;
}

inline int GrStat::getNrOfMacDeact() const
{
	return _nMacDeact;
}

#endif /* GRSTAT_H */
