/*
 * grsimulation.h
 *
 *  Created on: 05-nov-2009
 *      Author: M. El-Kebir
 */

#ifndef GRSIMULATION_H
#define GRSIMULATION_H

#include "gr.h"
#include "grsimulationgrid.h"
#include "grstat.h"
#include "macrophage.h"
#include "tgamma.h"
#include "tcytotoxic.h"
#include "tregulatory.h"
#include "ttest.h"
#include "grdiffusion.h"
#include "recruitmentbase.h"

class GrSimulation
{
private:
	int _time;
	GrSimulationGrid _grid;
	MacList _macList;
	TgamList _tgamList;
	TcytList _tcytList;
	TregList _tregList;
	GrStat _stats;
	double _areaThreshold;
	GrDiffusion* _pDiffusion;
	TTest* _pTTest[NOUTCOMES];
	RecruitmentBase* _pRecruitment;
	bool _tnfrDynamics;

	// Inhibits tnf secretion if true and if not using tnfr dynamics.
	bool _tnfKnockout;

	void moveTcells();
	void moveMacrophages();
	void updateStates();
	void updateT_Test();
	void computeNextStates();
	void secreteFromMacrophages();
	void secreteFromTcells();
	void secreteFromCaseations();
	void updateReceptorDynamics(double dt);
	void adjustTNFDegradation(double dt);
	void growExtMtb();
	void shuffleCells();

public:
	GrSimulation();
	~GrSimulation();
	void init();
	void solve();
	void performT_Test();
	int getTime() const;
	const GrStat& getStats() const;
	GrStat& getStats();
	const GrGrid& getGrid() const;
	GrGrid& getGrid();
	const MacList& getMacList() const;
	const TgamList& getTgamList() const;
	const TcytList& getTcytList() const;
	const TregList& getTregList() const;
	DiffusionMethod getDiffusionMethod() const;
	void setDiffusionMethod(DiffusionMethod method);
	bool getTnfrDynamics() const;
	void setTnfrDynamics(bool tnfrDynamics);
	double getAreaThreshold() const;
	void setAreaThreshold(double areaThreshold);
	OutcomeMethod getOutcomeMethod(int index) const;
	void getOutcomeParameters(int index, int& samplePeriod, int& testPeriod, double& alpha) const;
	void setOutcomeMethod(int index, OutcomeMethod method, double alpha, int testPeriod, int samplePeriod);
	void serialize(std::ostream& out) const;
	void deserialize(std::istream& in);
	void setRecruitment	(RecruitmentBase* pRecruitment);
	Mac* createMac(int row, int col, int birthtime, MacState state, bool NFkB, bool stat1);
	Tgam* createTgam(int row, int col, int birthtime, TgamState state);
	Tcyt* createTcyt(int row, int col, int birthtime, TcytState state);
	Treg* createTreg(int row, int col, int birthtime, TregState state);

	static void convertSimTime(const int time, int& rDays, int& rHours, int& rMinutes);
};

inline void GrSimulation::setRecruitment(RecruitmentBase* pRecruitment)
{
	_pRecruitment = pRecruitment;
}

inline void GrSimulation::getOutcomeParameters(int index,
	int& samplePeriod, int& testPeriod, double& alpha) const
{
	assert(0 <= index && index < NOUTCOMES);

	if (_pTTest[index])
	{
		samplePeriod = _pTTest[index]->getSamplePeriod();
		testPeriod = _pTTest[index]->getTestPeriod();
		alpha = _pTTest[index]->getAlpha();
	}
	else
	{
		samplePeriod = 1;
		testPeriod = 2;
		alpha = 0.05;
	}
}

inline OutcomeMethod GrSimulation::getOutcomeMethod(int index) const
{
	assert(0 <= index && index < NOUTCOMES);
	return _pTTest[index] ? _pTTest[index]->getMethod() : OUTCOME_NONE;
}

inline double GrSimulation::getAreaThreshold() const
{
	return _areaThreshold;
}

inline void GrSimulation::setAreaThreshold(double areaThreshold)
{
	_areaThreshold = areaThreshold;
}

inline void GrSimulation::convertSimTime(const int time, int& rDays, int& rHours, int& rMinutes)
{
	rDays = time / 144;
	rHours = (time % 144) / 6;
	rMinutes = ((time % 144) % 6) * 10;
}

inline int GrSimulation::getTime() const
{
	return _time;
}

inline DiffusionMethod GrSimulation::getDiffusionMethod() const
{
	return _pDiffusion->getMethod();
}

inline bool GrSimulation::getTnfrDynamics() const
{
	return _tnfrDynamics;
}

inline void GrSimulation::setTnfrDynamics(bool tnfrDynamics)
{
	_tnfrDynamics = tnfrDynamics;
}

inline const MacList& GrSimulation::getMacList() const
{
	return _macList;
}

inline const TgamList& GrSimulation::getTgamList() const
{
	return _tgamList;
}

inline const TcytList& GrSimulation::getTcytList() const
{
	return _tcytList;
}

inline const TregList& GrSimulation::getTregList() const
{
	return _tregList;
}

inline const GrGrid& GrSimulation::getGrid() const
{
	return _grid.getGrid();
}

inline GrGrid& GrSimulation::getGrid()
{
	return _grid.getGrid();
}

inline const GrStat& GrSimulation::getStats() const
{
	return _stats;
}

inline GrStat& GrSimulation::getStats()
{
	return _stats;
}

inline Mac* GrSimulation::createMac(int row, int col, int birthtime, MacState state, bool NFkB, bool stat1)
{
	_macList.push_back(Mac(birthtime, row, col, state, 0, NFkB, stat1));
	Mac* pMac = &_macList.back();
	
	assert_res(_grid(row, col).addAgent(pMac));
	_stats.updateMacStatistics(state);

	return pMac;
}

inline Tgam* GrSimulation::createTgam(int row, int col, int birthtime, TgamState state)
{
	_tgamList.push_back(Tgam(birthtime, row, col, state));
	Tgam* pTgam = &_tgamList.back();
	
	assert_res(_grid(row, col).addAgent(pTgam));
	_stats.updateTgamStatistics(state);

	return pTgam;
}

inline Tcyt* GrSimulation::createTcyt(int row, int col, int birthtime, TcytState state)
{
	_tcytList.push_back(Tcyt(birthtime, row, col, state));
	Tcyt* pTcyt = &_tcytList.back();
	
	assert_res(_grid(row, col).addAgent(pTcyt));
	_stats.updateTcytStatistics(state);

	return pTcyt;
}

inline Treg* GrSimulation::createTreg(int row, int col, int birthtime, TregState state)
{
	_tregList.push_back(Treg(birthtime, row, col, state));
	Treg* pTreg = &_tregList.back();
	
	assert_res(_grid(row, col).addAgent(pTreg));
	_stats.updateTregStatistics(state);

	return pTreg;
}

#endif /* GRSIMULATION_H */
