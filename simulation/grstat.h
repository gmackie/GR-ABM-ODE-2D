/*
 * grstat.h
 *
 *  Created on: 30-nov-2009
 *      Author: M. El-Kebir
 */

#ifndef GRSTAT_H
#define GRSTAT_H

#include "gr.h"
#include "params.h"
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/median.hpp>

typedef enum {	GR_CONTAINMENT, GR_CONTAINMENT_INCONSISTENT, GR_CLEARANCE,
				GR_DISSEMINATION, GR_DISSEMINATION_INCONSISTENT, GR_UNKNOWN, GR_NONE} GrStatus;

namespace ba = boost::accumulators;
class GrStat
{
private:
	static const std::string _ClassName;

	/*
	 * !!! If the data members change then the serialize and deserialize functions need to be updated !!!
	 */

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

	double _totExtMtb;			// Replicating and non-replicating bacteria: all ExtMtb in all compartments
	double _totNonRepExtMtb;	// Non-replicating bacteria: in caseated compartments.
	double _totIntMtb;

	double _totMacAttractant;
	double _totTNF;
	double _totCCL2;
	double _totCCL5;
	double _totCXCL9;

	int _nApoptosisFasFasL;
	int _nMacApoptosisTNF;
	int _nRestingMacApoptosisTNF;
	int _nInfAndCinfMacApoptosisTNF;
	int _nActivatedMacApoptosisTNF;
	int _nTcellApoptosisTNF;
	
	int _nRestingMacActivationTNF;
	int _nInfMacActivationTNF;

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
	int _areaCellDensity;
	int _nCaseated;

	GrStatus _grStatus[NOUTCOMES];

	int _queueTgam;
	int _queueTcyt;
	int _queueTreg;
	double _fluxTgam;
	double _fluxTcyt;
	double _fluxTreg;
	int _nSourceMacActive;
	int _nSourceTgamActive;
	int _nSourceTcytActive;
	int _nSourceTregActive;
	double _MDC;
	double _N4;
	double _TH0;
	double _TH1;
	double _N8;
	double _T80;
	double _T8;
	double _TC;
	double _TH0lung;
	double _TH1lung;
	double _T80lung;
	double _T8lung;
	double _TClung;

  unsigned* _intMtbFreq;
public:
  typedef ba::accumulator_set<double, ba::stats< ba::features< ba::tag::variance, ba::tag::min, ba::tag::max, ba::tag::median > > > Stat;
private:
  Stat* _macIntMtbStats;

public:
	GrStat();
	~GrStat();
  const unsigned* getIntMtbFreq(size_t& s) const;
  const Stat* getIntMtbStats(size_t &s) const;
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
	double getTotNonRepExtMtb() const;
	double getTotIntMtb() const;
	double getTotMacAttractant() const;
	double getTotTNF() const;
	double getTotCCL2() const;
	double getTotCCL5() const;
	double getTotCXCL9() const;
	void updateAgentStatistics(Agent* a);
protected:
	void updateMacStatistics(MacState state);
	void updateTgamStatistics(TgamState state);
	void updateTcytStatistics(TcytState state);
	void updateTregStatistics(TregState state);
public:
	void updateMacNFkBStatistics(MacState state);
	void updateMacStat1Statistics(MacState state);
	void updateMacDeactStatistics(MacState state);
	void resetAgentStats();
	void reset();
	void incTotExtMtb(double dExtMtb);
	void incTotNonRepExtMtb(double dNonRepExtMtb);
	void incTotIntMtb(double dIntMtb);
	void incTotMacAttractant(double dMacAttractant);
	void incTotTNF(double dTNF);
	void incTotCCL2(double dCCL2);
	void incTotCCL5(double dCCL5);
	void incTotCXCL9(double dCXCL9);
	int getNrApoptosisFasFasL() const;
	int getNrMacApoptosisTNF() const;
	int getNrRestingMacApoptosisTNF() const;
	int getNrInfAndCinfMacApoptosisTNF() const;
	int getNrActivatedMacApoptosisTNF() const;
	int getNrTcellApoptosisTNF() const;
	int getNrRestingMacActivationTNF() const;
	int getNrInfMacActivationTNF() const;
	void incApoptosisFasFasL();
	void incMacApoptosisTNF();
	void incRestingMacApoptosisTNF();
	void incInfAndCinfMacApoptosisTNF();
	void incActivatedMacApoptosisTNF();
	void incTcellApoptosisTNF();
	void incRestingMacActivationTNF();
	void incInfMacActivationTNF();
	int getNrSourcesMac() const;
	int getNrSourcesTgam() const;
	int getNrSourcesTcyt() const;
	int getNrSourcesTreg() const;
	void incNrSourcesMac();
	void incNrSourcesTgam();
	void incNrSourcesTcyt();
	void incNrSourcesTreg();
	int getNrSourcesActiveMac() const;
	int getNrSourcesActiveTgam() const;
	int getNrSourcesActiveTcyt() const;
	int getNrSourcesActiveTreg() const;
	void incNrSourcesActiveMac();
	void incNrSourcesActiveTgam();
	void incNrSourcesActiveTcyt();
	void incNrSourcesActiveTreg();
	int getNrBactAct() const;
	void incNrBactAct();
	int getArea() const;
	void incArea();
	int getAreaCellDensity() const;
	void incAreaCellDensity();
	int getNrCaseated() const;
	void incNrCaseated();
	GrStatus getGrStatus(int index) const;
	void setGrStatus(int index, GrStatus status);
	void serialize(std::ostream& out) const;
	void deserialize(std::istream& in);
	int getNrTgamQueued() const;
	void setNrTgamQueued(int count);
	int getNrTcytQueued() const;
	void setNrTcytQueued(int count);
	int getNrTregQueued() const;
	void setNrTregQueued(int count);
	double getFluxTgam() const;
	void setFluxTgam(double flux);
	double getFluxTcyt() const;
	void setFluxTcyt(double flux);
	double getFluxTreg() const;
	void setFluxTreg(double flux);
	double getMDC() const;
	void setMDC(double MDC);
	double getN4() const;
	double getTH0() const;
	double getTH1() const;
	double getN8() const;
	double getT80() const;
	double getT8() const;
	double getTC() const;
	double getTH0lung() const;
	double getTH1lung() const;
	double getT80lung() const;
	double getT8lung() const;
	double getTClung() const;
	void setN4(double val);
	void setTH0(double val);
	void setTH1(double val);
	void setN8(double val);
	void setT80(double val);
	void setT8(double val);
	void setTC(double val);
	void setTH0lung(double val);
	void setTH1lung(double val);
	void setT80lung(double val);
	void setT8lung(double val);
	void setTClung(double val);
};

inline const unsigned* GrStat::getIntMtbFreq(size_t& s) const {
  s = size_t(_PARAM(PARAM_MAC_THRESHOLD_BURST_CI_INTMTB));
  return _intMtbFreq;
}
inline const GrStat::Stat* GrStat::getIntMtbStats(size_t &s) const {
  s = size_t(NMAC_STATES);
  return _macIntMtbStats;
}
inline void GrStat::setN4(double val)
{
	_N4 = val;
}

inline void GrStat::setTH0(double val)
{
	_TH0 = val;
}

inline void GrStat::setTH1(double val)
{
	_TH1 = val;
}

inline void GrStat::setN8(double val)
{
	_N8 = val;
}

inline void GrStat::setT80(double val)
{
	_T80 = val;
}

inline void GrStat::setT8(double val)
{
	_T8 = val;
}

inline void GrStat::setTC(double val)
{
	_TC = val;
}

inline void GrStat::setTH0lung(double val)
{
	_TH0lung = val;
}

inline void GrStat::setTH1lung(double val)
{
	_TH1lung = val;
}

inline void GrStat::setT80lung(double val)
{
	_T80lung = val;
}

inline void GrStat::setT8lung(double val)
{
	_T8lung = val;
}

inline void GrStat::setTClung(double val)
{
	_TClung = val;
}

inline double GrStat::getN4() const
{
	return _N4;
}

inline double GrStat::getTH0() const
{
	return _TH0;
}

inline double GrStat::getTH1() const
{
	return _TH1;
}

inline double GrStat::getN8() const
{
	return _N8;
}

inline double GrStat::getT80() const
{
	return _T80;
}

inline double GrStat::getT8() const
{
	return _T8;
}

inline double GrStat::getTC() const
{
	return _TC;
}

inline double GrStat::getTH0lung() const
{
	return _TH0lung;
}

inline double GrStat::getTH1lung() const
{
	return _TH1lung;
}

inline double GrStat::getT80lung() const
{
	return _T80lung;
}

inline double GrStat::getT8lung() const
{
	return _T8lung;
}

inline double GrStat::getTClung() const
{
	return _TClung;
}

inline double GrStat::getMDC() const
{
	return _MDC;
}

inline void GrStat::setMDC(double MDC)
{
	_MDC = MDC;
}

inline int GrStat::getNrTgamQueued() const
{
	return _queueTgam;
}

inline void GrStat::setNrTgamQueued(int count)
{
	_queueTgam = count;
}

inline double GrStat::getFluxTgam() const
{
	return _fluxTgam;
}

inline void GrStat::setFluxTgam(double flux)
{
	_fluxTgam = flux;
}

inline int GrStat::getNrTcytQueued() const
{
	return _queueTcyt;
}

inline void GrStat::setNrTcytQueued(int count)
{
	_queueTcyt = count;
}

inline double GrStat::getFluxTcyt() const
{
	return _fluxTcyt;
}

inline void GrStat::setFluxTcyt(double flux)
{
	_fluxTcyt = flux;
}

inline int GrStat::getNrTregQueued() const
{
	return _queueTreg;
}

inline void GrStat::setNrTregQueued(int count)
{
	_queueTreg = count;
}

inline double GrStat::getFluxTreg() const
{
	return _fluxTreg;
}

inline void GrStat::setFluxTreg(double flux)
{
	_fluxTreg = flux;
}

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

inline int GrStat::getAreaCellDensity() const
{
	return _areaCellDensity;
}

inline void GrStat::incAreaCellDensity()
{
	_areaCellDensity++;
}

inline int GrStat::getNrCaseated() const
{
	return _nCaseated;
}

inline void GrStat::incNrCaseated()
{
	_nCaseated++;
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

inline int GrStat::getNrSourcesActiveMac() const
{
	return _nSourceMacActive;
}

inline int GrStat::getNrSourcesActiveTgam() const
{
	return _nSourceTgamActive;
}

inline int GrStat::getNrSourcesActiveTcyt() const
{
	return _nSourceTcytActive;
}

inline int GrStat::getNrSourcesActiveTreg() const
{
	return _nSourceTregActive;
}

inline void GrStat::incNrSourcesActiveMac()
{
	_nSourceMacActive++;
}

inline void GrStat::incNrSourcesActiveTgam()
{
	_nSourceTgamActive++;
}

inline void GrStat::incNrSourcesActiveTcyt()
{
	_nSourceTcytActive++;
}

inline void GrStat::incNrSourcesActiveTreg()
{
	_nSourceTregActive++;
}

inline int GrStat::getNrApoptosisFasFasL() const
{
	return _nApoptosisFasFasL;
}

inline int GrStat::getNrMacApoptosisTNF() const
{
	return _nMacApoptosisTNF;
}

inline int GrStat::getNrRestingMacApoptosisTNF() const
{
	return _nRestingMacApoptosisTNF;
}

inline int GrStat::getNrInfAndCinfMacApoptosisTNF() const
{
	return _nInfAndCinfMacApoptosisTNF;
}

inline int GrStat::getNrActivatedMacApoptosisTNF() const
{
	return _nActivatedMacApoptosisTNF;
}

inline int GrStat::getNrTcellApoptosisTNF() const
{
	return _nTcellApoptosisTNF;
}

inline int GrStat::getNrRestingMacActivationTNF() const
{
	return _nRestingMacActivationTNF;
}

inline int GrStat::getNrInfMacActivationTNF() const
{
	return _nInfMacActivationTNF;
}

inline void GrStat::incApoptosisFasFasL()
{
	_nApoptosisFasFasL++;
}

inline void GrStat::incMacApoptosisTNF()
{
	_nMacApoptosisTNF++;
}

inline void GrStat::incRestingMacApoptosisTNF()
{
	_nRestingMacApoptosisTNF++;
}

inline void GrStat::incInfAndCinfMacApoptosisTNF()
{
	_nInfAndCinfMacApoptosisTNF++;
}

inline void GrStat::incActivatedMacApoptosisTNF()
{
	_nActivatedMacApoptosisTNF++;
}

inline void GrStat::incTcellApoptosisTNF()
{
	_nTcellApoptosisTNF++;
}

inline void GrStat::incRestingMacActivationTNF()
{
	_nRestingMacActivationTNF++;
}

inline void GrStat::incInfMacActivationTNF()
{
	_nInfMacActivationTNF++;
}

inline void GrStat::incTotExtMtb(double dExtMtb)
{
	_totExtMtb += dExtMtb;
}

inline void GrStat::incTotNonRepExtMtb(double dNonRepExtMtb)
{
	_totNonRepExtMtb += dNonRepExtMtb;
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

inline double GrStat::getTotNonRepExtMtb() const
{
	return _totNonRepExtMtb;
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
