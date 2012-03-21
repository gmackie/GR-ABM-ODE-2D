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
#include <functional>
#include <algorithm>
#include <numeric>
#include <sstream>


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

	int _nAgents[NAGENTS];
	int _nMac[NMAC_STATES];
	int _nTgam[NTGAM_STATES];
	int _nTcyt[NTCYT_STATES];
	int _nTreg[NTREG_STATES];

	Scalar _totExtMtb;			// Replicating and non-replicating bacteria: all ExtMtb in all compartments
	Scalar _totNonRepExtMtb;	// Non-replicating bacteria: in caseated compartments.
	Scalar _totIntMtb;

	Scalar _totMacAttractant;
	Scalar _totTNF;
    Scalar _totIL10;
	Scalar _totCCL2;
	Scalar _totCCL5;
	Scalar _totCXCL9;

    Scalar _totTNFR1int;
    Scalar _totkmRNA;

	int _nApoptosisFasFasL;
	int _nMacApoptosisTNF[NMAC_STATES];
	int _nTcellApoptosisTNF;
	
	int _nRestingMacActivationTNF;
	int _nInfMacActivationTNF;

	int _nMacNFkB[NMAC_STATES];

	int _nSource[NAGENTS];
	int _nSourceActive[NAGENTS];

	int _nSourceMacCrowded;	// Should be an array
	int _nSourceTgamCrowded;
	int _nSourceTcytCrowded;
	int _nSourceTregCrowded;

	int _nMacStat1[NMAC_STATES];
	int _nMacDeact[NMAC_STATES];

	int _nBactAct;
	int _areaTNF;
	int _areaCellDensity;
	int _nCaseated;

	GrStatus _grStatus[NOUTCOMES];

	int _queueTgam;
	int _queueTcyt;
	int _queueTreg;

	int _queueTgamDie;
	int _queueTcytDie;
	int _queueTregDie;

	int _recruitedTgam;
	int _recruitedTcyt;
	int _recruitedTreg;

	Scalar _fluxTgam;
	Scalar _fluxTcyt;
	Scalar _fluxTreg;

	Scalar _MDC;
	Scalar _N4;
	Scalar _TH0;
	Scalar _TH1;
	Scalar _N8;
	Scalar _T80;
	Scalar _T8;
	Scalar _TC;
	Scalar _TH0lung;
	Scalar _TH1lung;
	Scalar _T80lung;
	Scalar _T8lung;
	Scalar _TClung;
	int _nCellTnfInhibit;
	
	size_t _intMtbFreqSize;
	unsigned* _intMtbFreq;
public:
  typedef ba::accumulator_set<Scalar, ba::stats< ba::features< ba::tag::variance, ba::tag::min, ba::tag::max, ba::tag::median > > > Stat;
private:
  Stat _macIntMtbStats[NMAC_STATES];

public:
	GrStat();
  GrStat(const GrStat&);
  GrStat& operator=(const GrStat&);
	~GrStat();
	size_t getIntMtbFreqSize() const;
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
    int getNrOfTgamDouble() const;
    int getNrOfTgamInduced() const;
	int getNrOfTcytActive() const;
	int getNrOfTcytDead() const;
	int getNrOfTcytDownRegulated() const;
	int getNrOfTregActive() const;
	int getNrOfTregDead() const;
	Scalar getTotExtMtb() const;
	Scalar getTotNonRepExtMtb() const;
	Scalar getTotIntMtb() const;
	Scalar getTotMacAttractant() const;
	Scalar getTotTNF() const;
    Scalar getTotIL10() const;
	Scalar getTotCCL2() const;
	Scalar getTotCCL5() const;
	Scalar getTotCXCL9() const;
    Scalar getTotTNFR1int() const;
    Scalar getTotkmRNA() const;
	void updateAgentStatistics(Agent* a);
public:
	void updateMacNFkBStatistics(MacState state);
	void updateMacStat1Statistics(MacState state);
	void updateMacDeactStatistics(MacState state);
	void resetAgentStats();
	void reset();
	void incTotExtMtb(Scalar dExtMtb);
	void incTotNonRepExtMtb(Scalar dNonRepExtMtb);
	void incTotIntMtb(Scalar dIntMtb);
	void incTotMacAttractant(Scalar dMacAttractant);
	void incTotTNF(Scalar dTNF);
    void incTotIL10(Scalar dIL10);
	void incTotCCL2(Scalar dCCL2);
	void incTotCCL5(Scalar dCCL5);
	void incTotCXCL9(Scalar dCXCL9);
    void incTotTNFR1int(Scalar dTNFR1int);
    void incTotkmRNA(Scalar dkmRNA);
	int getNrApoptosisFasFasL() const;
	int getNrMacApoptosisTNF() const;
	int getNrMacApoptosisTNF(MacState s) const;
	int getNrTcellApoptosisTNF() const;
	int getNrRestingMacActivationTNF() const;
	int getNrInfMacActivationTNF() const;
	void incApoptosisFasFasL();
	void incMacApoptosisTNF(MacState s);
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

	int getNrSourcesCrowdedMac() const;
	int getNrSourcesCrowdedTgam() const;
	int getNrSourcesCrowdedTcyt() const;
	int getNrSourcesCrowdedTreg() const;
	void incNrSourcesCrowdedMac();
	void incNrSourcesCrowdedTgam();
	void incNrSourcesCrowdedTcyt();
	void incNrSourcesCrowdedTreg();


	int getNrBactAct() const;
	void incNrBactAct();
	int getAreaTNF() const;
	void incAreaTNF();
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
	int getNrTgamQueuedDie() const;
	void incNrTgamQueuedDie();
	int getNrTcytQueuedDie() const;
	void incNrTcytQueuedDie();
	int getNrTregQueuedDie() const;
	void incNrTregQueuedDie();
	void incNrQueuedDie(TcellType type);

	int getNrTgamRecruited() const;
	void incNrTgamRecruited();
	int getNrTcytRecruited() const;
	void incNrTcytRecruited();
	int getNrTregRecruited() const;
	void incNrTregRecruited();

	Scalar getFluxTgam() const;
	void setFluxTgam(Scalar flux);
	Scalar getFluxTcyt() const;
	void setFluxTcyt(Scalar flux);
	Scalar getFluxTreg() const;
	void setFluxTreg(Scalar flux);
	Scalar getMDC() const;
	void setMDC(Scalar MDC);
	Scalar getN4() const;
	Scalar getTH0() const;
	Scalar getTH1() const;
	Scalar getN8() const;
	Scalar getT80() const;
	Scalar getT8() const;
	Scalar getTC() const;
	Scalar getTH0lung() const;
	Scalar getTH1lung() const;
	Scalar getT80lung() const;
	Scalar getT8lung() const;
	Scalar getTClung() const;
	void setN4(Scalar val);
	void setTH0(Scalar val);
	void setTH1(Scalar val);
	void setN8(Scalar val);
	void setT80(Scalar val);
	void setT8(Scalar val);
	void setTC(Scalar val);
	void setTH0lung(Scalar val);
	void setTH1lung(Scalar val);
	void setT80lung(Scalar val);
	void setT8lung(Scalar val);
	void setTClung(Scalar val);
	int getNrOfCellsInhibited () const;
	void incNrOfCellsInhibited ();
};

inline size_t GrStat::getIntMtbFreqSize() const
{
	return _intMtbFreqSize;
}

inline const unsigned* GrStat::getIntMtbFreq(size_t& s) const {
  s = size_t(_PARAM(PARAM_MAC_THRESHOLD_BURST_CI_INTMTB));
  return _intMtbFreq;
}
inline const GrStat::Stat* GrStat::getIntMtbStats(size_t &s) const {
  s = size_t(NMAC_STATES);
  return _macIntMtbStats;
}
inline void GrStat::setN4(Scalar val)
{
	_N4 = val;
}

inline void GrStat::setTH0(Scalar val)
{
	_TH0 = val;
}

inline void GrStat::setTH1(Scalar val)
{
	_TH1 = val;
}

inline void GrStat::setN8(Scalar val)
{
	_N8 = val;
}

inline void GrStat::setT80(Scalar val)
{
	_T80 = val;
}

inline void GrStat::setT8(Scalar val)
{
	_T8 = val;
}

inline void GrStat::setTC(Scalar val)
{
	_TC = val;
}

inline void GrStat::setTH0lung(Scalar val)
{
	_TH0lung = val;
}

inline void GrStat::setTH1lung(Scalar val)
{
	_TH1lung = val;
}

inline void GrStat::setT80lung(Scalar val)
{
	_T80lung = val;
}

inline void GrStat::setT8lung(Scalar val)
{
	_T8lung = val;
}

inline void GrStat::setTClung(Scalar val)
{
	_TClung = val;
}

inline Scalar GrStat::getN4() const
{
	return _N4;
}

inline Scalar GrStat::getTH0() const
{
	return _TH0;
}

inline Scalar GrStat::getTH1() const
{
	return _TH1;
}

inline Scalar GrStat::getN8() const
{
	return _N8;
}

inline Scalar GrStat::getT80() const
{
	return _T80;
}

inline Scalar GrStat::getT8() const
{
	return _T8;
}

inline Scalar GrStat::getTC() const
{
	return _TC;
}

inline Scalar GrStat::getTH0lung() const
{
	return _TH0lung;
}

inline Scalar GrStat::getTH1lung() const
{
	return _TH1lung;
}

inline Scalar GrStat::getT80lung() const
{
	return _T80lung;
}

inline Scalar GrStat::getT8lung() const
{
	return _T8lung;
}

inline Scalar GrStat::getTClung() const
{
	return _TClung;
}

inline Scalar GrStat::getMDC() const
{
	return _MDC;
}

inline void GrStat::setMDC(Scalar MDC)
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

inline int GrStat::getNrTgamQueuedDie() const
{
	return _queueTgamDie;
}

inline void GrStat::incNrTgamQueuedDie()
{
	_queueTgamDie++;
}

inline int GrStat::getNrTgamRecruited() const
{
	return _recruitedTgam;
}

inline void GrStat::incNrTgamRecruited()
{
	_recruitedTgam++;
}

inline Scalar GrStat::getFluxTgam() const
{
	return _fluxTgam;
}

inline void GrStat::setFluxTgam(Scalar flux)
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

inline int GrStat::getNrTcytQueuedDie() const
{
	return _queueTcytDie;
}

inline void GrStat::incNrTcytQueuedDie()
{
	_queueTcytDie++;
}

inline int GrStat::getNrTcytRecruited() const
{
	return _recruitedTcyt;
}

inline void GrStat::incNrTcytRecruited()
{
	_recruitedTcyt++;
}

inline Scalar GrStat::getFluxTcyt() const
{
	return _fluxTcyt;
}

inline void GrStat::setFluxTcyt(Scalar flux)
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

inline int GrStat::getNrTregQueuedDie() const
{
	return _queueTregDie;
}

inline void GrStat::incNrTregQueuedDie()
{
	_queueTregDie++;
}

inline void GrStat::incNrQueuedDie(TcellType type)
{
	switch (type)
	{
		case TCELL_TYPE_CYT:
			incNrTcytQueuedDie();
			break;
		case TCELL_TYPE_REG:
			incNrTregQueuedDie();
			break;
		case TCELL_TYPE_GAM:
			incNrTgamQueuedDie();
			break;
		default:
			std::ostringstream os;
			os << "GrStat::incNrQueuedDie: invalid T cell type of " << type << ", cannot continue...";
			throw std::runtime_error(os.str());
	}
}

inline int GrStat::getNrTregRecruited() const
{
	return _recruitedTreg;
}

inline void GrStat::incNrTregRecruited()
{
	_recruitedTreg++;
}

inline Scalar GrStat::getFluxTreg() const
{
	return _fluxTreg;
}

inline void GrStat::setFluxTreg(Scalar flux)
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

inline int GrStat::getAreaTNF() const
{
	return _areaTNF;
}

inline void GrStat::incAreaTNF()
{
	_areaTNF++;
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
	return _nSource[MAC];
}

inline int GrStat::getNrSourcesTgam() const
{
	return _nSource[TGAM];
}

inline int GrStat::getNrSourcesTcyt() const
{
	return _nSource[TCYT];
}

inline int GrStat::getNrSourcesTreg() const
{
	return _nSource[TREG];
}

inline void GrStat::incNrSourcesMac()
{
	_nSource[MAC]++;
}

inline void GrStat::incNrSourcesTgam()
{
	_nSource[TGAM]++;
}

inline void GrStat::incNrSourcesTcyt()
{
	_nSource[TCYT]++;
}

inline void GrStat::incNrSourcesTreg()
{
	_nSource[TREG]++;
}

inline int GrStat::getNrSourcesActiveMac() const
{
	return _nSourceActive[MAC];
}

inline int GrStat::getNrSourcesActiveTgam() const
{
	return _nSourceActive[TGAM];
}

inline int GrStat::getNrSourcesActiveTcyt() const
{
	return _nSourceActive[TCYT];
}

inline int GrStat::getNrSourcesActiveTreg() const
{
	return _nSourceActive[TREG];
}

inline void GrStat::incNrSourcesActiveMac()
{
	_nSourceActive[MAC]++;
}

inline void GrStat::incNrSourcesActiveTgam()
{
	_nSourceActive[TGAM]++;
}

inline void GrStat::incNrSourcesActiveTcyt()
{
	_nSourceActive[TCYT]++;
}

inline void GrStat::incNrSourcesActiveTreg()
{
	_nSourceActive[TREG]++;
}

inline int GrStat::getNrSourcesCrowdedMac() const
{
	return _nSourceMacCrowded;
}

inline int GrStat::getNrSourcesCrowdedTgam() const
{
	return _nSourceTgamCrowded;
}

inline int GrStat::getNrSourcesCrowdedTcyt() const
{
	return _nSourceTcytCrowded;
}

inline int GrStat::getNrSourcesCrowdedTreg() const
{
	return _nSourceTregCrowded;
}

inline void GrStat::incNrSourcesCrowdedMac()
{
	_nSourceMacCrowded++;
}

inline void GrStat::incNrSourcesCrowdedTgam()
{
	_nSourceTgamCrowded++;
}

inline void GrStat::incNrSourcesCrowdedTcyt()
{
	_nSourceTcytCrowded++;
}

inline void GrStat::incNrSourcesCrowdedTreg()
{
	_nSourceTregCrowded++;
}

inline int GrStat::getNrApoptosisFasFasL() const
{
	return _nApoptosisFasFasL;
}

inline int GrStat::getNrMacApoptosisTNF() const
{
  return std::accumulate(_nMacApoptosisTNF, _nMacApoptosisTNF+NMAC_STATES, 0);
}
inline int GrStat::getNrMacApoptosisTNF(MacState s) const
{
  return _nMacApoptosisTNF[s];
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

inline void GrStat::incMacApoptosisTNF(MacState s)
{
	_nMacApoptosisTNF[s]++;
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

inline void GrStat::incTotExtMtb(Scalar dExtMtb)
{
	_totExtMtb += dExtMtb;
}

inline void GrStat::incTotNonRepExtMtb(Scalar dNonRepExtMtb)
{
	_totNonRepExtMtb += dNonRepExtMtb;
}

inline void GrStat::incTotIntMtb(Scalar dIntMtb)
{
	_totIntMtb += dIntMtb;
}

inline void GrStat::incTotMacAttractant(Scalar dMacAttractant)
{
	_totMacAttractant += dMacAttractant;
}

inline void GrStat::incTotTNF(Scalar dTNF)
{
	_totTNF += dTNF;
}

inline void GrStat::incTotIL10(Scalar dIL10)
{
    _totIL10 += dIL10;
}

inline void GrStat::incTotCCL2(Scalar dCCL2)
{
	_totCCL2 += dCCL2;
}

inline void GrStat::incTotCCL5(Scalar dCCL5)
{
	_totCCL5 += dCCL5;
}

inline void GrStat::incTotCXCL9(Scalar dCXCL9)
{
	_totCXCL9 += dCXCL9;
}

inline void GrStat::incTotTNFR1int(Scalar dTNFR1int)
{
    _totTNFR1int += dTNFR1int;
}

inline void GrStat::incTotkmRNA(Scalar dkmRNA)
{
    _totkmRNA += dkmRNA;
}

inline Scalar GrStat::getTotMacAttractant() const
{
	return _totMacAttractant;
}

inline Scalar GrStat::getTotTNF() const
{
	return _totTNF;
}

inline Scalar GrStat::getTotIL10() const
{
    return _totIL10;
}

inline Scalar GrStat::getTotCCL2() const
{
	return _totCCL2;
}

inline Scalar GrStat::getTotCCL5() const
{
	return _totCCL5;
}

inline Scalar GrStat::getTotCXCL9() const
{
	return _totCXCL9;
}

inline Scalar GrStat::getTotTNFR1int() const
{
    return _totTNFR1int;
}

inline Scalar GrStat::getTotkmRNA() const
{
    return _totkmRNA;
}

inline Scalar GrStat::getTotExtMtb() const
{
	return _totExtMtb;
}

inline Scalar GrStat::getTotNonRepExtMtb() const
{
	return _totNonRepExtMtb;
}

inline Scalar GrStat::getTotIntMtb() const
{
	return _totIntMtb;
}

inline int GrStat::getNrOfMacResting() const
{
	return _nMac[MAC_RESTING];
}

inline int GrStat::getNrOfMacInfected() const
{
	return _nMac[MAC_INFECTED];
}

inline int GrStat::getNrOfMacCInfected() const
{
	return _nMac[MAC_CINFECTED];
}

inline int GrStat::getNrOfMacActive() const
{
	return _nMac[MAC_ACTIVE];
}

inline int GrStat::getNrOfMacDead() const
{
	return _nMac[MAC_DEAD];
}

inline int GrStat::getNrOfMac() const
{
	return _nAgents[MAC];
}

inline int GrStat::getNrOfMacNFkBResting() const
{
	return _nMacNFkB[MAC_RESTING];
}

inline int GrStat::getNrOfMacNFkBInfected() const
{
	return _nMacNFkB[MAC_INFECTED];
}

inline int GrStat::getNrOfMacNFkBCInfected() const
{
	return _nMacNFkB[MAC_CINFECTED];
}

inline int GrStat::getNrOfMacNFkBActive() const
{
	return _nMacNFkB[MAC_ACTIVE];
}

inline int GrStat::getNrOfMacNFkBDead() const
{
	return _nMacNFkB[MAC_DEAD];
}

inline int GrStat::getNrOfMacNFkB() const
{
	return std::accumulate(_nMacNFkB, _nMacNFkB+NMAC_STATES, 0);
}

inline int GrStat::getNrOfTgamActive() const
{
	return _nTgam[TGAM_ACTIVE];
}

inline int GrStat::getNrOfTgamDead() const
{
	return _nTgam[TGAM_DEAD];
}

inline int GrStat::getNrOfTgamDownRegulated() const
{
	return _nTgam[TGAM_DOWN_REGULATED];
}

inline int GrStat::getNrOfTgam() const
{
  return _nAgents[TGAM];
}

inline int GrStat::getNrOfTgamDouble() const
{
    return _nTgam[TGAM_ACTIVE_DOUBLE];
}

inline int GrStat::getNrOfTgamInduced() const
{
    return _nTgam[TGAM_INDUCED_REG];
}

inline int GrStat::getNrOfTcytDead() const
{
	return _nTcyt[TCYT_DEAD];
}

inline int GrStat::getNrOfTcytActive() const
{
	return _nTcyt[TCYT_ACTIVE];
}

inline int GrStat::getNrOfTcytDownRegulated() const
{
	return _nTcyt[TCYT_DOWN_REGULATED];
}

inline int GrStat::getNrOfTcyt() const
{
	return _nAgents[TCYT];
}

inline int GrStat::getNrOfTregDead() const
{
	return _nTreg[TREG_DEAD];
}

inline int GrStat::getNrOfTregActive() const
{
	return _nTreg[TREG_ACTIVE];
}

inline int GrStat::getNrOfTreg() const
{
	return _nAgents[TREG];
}

inline int GrStat::getNrOfMacStat1Resting() const
{
	return _nMacStat1[MAC_RESTING];
}

inline int GrStat::getNrOfMacStat1Infected() const
{
	return _nMacStat1[MAC_INFECTED];
}

inline int GrStat::getNrOfMacStat1CInfected() const
{
	return _nMacStat1[MAC_CINFECTED];
}

inline int GrStat::getNrOfMacStat1Active() const
{
	return _nMacStat1[MAC_ACTIVE];
}

inline int GrStat::getNrOfMacStat1Dead() const
{
	return _nMacStat1[MAC_DEAD];
}

inline int GrStat::getNrOfMacStat1() const
{
	return std::accumulate(_nMacStat1,_nMacStat1+NMAC_STATES,0);
}

inline int GrStat::getNrOfMacDeactResting() const
{
	return _nMacDeact[MAC_RESTING];
}

inline int GrStat::getNrOfMacDeactInfected() const
{
	return _nMacDeact[MAC_INFECTED];
}

inline int GrStat::getNrOfMacDeactCInfected() const
{
	return _nMacDeact[MAC_CINFECTED];
}

inline int GrStat::getNrOfMacDeactActive() const
{
	return _nMacDeact[MAC_ACTIVE];
}

inline int GrStat::getNrOfMacDeactDead() const
{
	return _nMacDeact[MAC_DEAD];
}

inline int GrStat::getNrOfMacDeact() const
{
	return std::accumulate(_nMacDeact,_nMacDeact+NMAC_STATES,0);
}

inline int GrStat::getNrOfCellsInhibited() const
{
	return _nCellTnfInhibit;
}

inline void GrStat::incNrOfCellsInhibited()
{
	_nCellTnfInhibit++;
}


#endif /* GRSTAT_H */
