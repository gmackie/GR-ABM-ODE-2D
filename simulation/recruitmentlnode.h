/*
 * recruitmentlnode.h
 *
 *  Created on: Jun 17, 2010
 *      Author: mohammed
 */

#ifndef RECRUITMENTLNODE_H_
#define RECRUITMENTLNODE_H_

#include "gr.h"
#include "grgrid.h"
#include "grstat.h"
#include "recruitmentbase.h"
#include <string>

class RecruitmentLnODE : public RecruitmentBase
{
protected:
	enum TcellType {TCELL_TYPE_CYT, TCELL_TYPE_REG, TCELL_TYPE_GAM, TCELL_TYPE_COUNT};
	static const int _nrConditions = 43;

	static const int _idxMDC = 13;
	static const int _idxNaiveCD4 = 14;
	static const int _idxNaiveCD8 = 16;
	static const int _idxEffectorTH1 = 39;
	static const int _idxEffectorT8 = 41;
	static const int _idxCTL = 42;

	double _tcellTable[TCELL_TYPE_COUNT]; // contains lower bounds

	int _tcellQueueCount[TCELL_TYPE_COUNT];
	std::vector<TcellType> _tcellQueue;

	int _prevMiMci;
	double _odeInitialConditions[_nrConditions];
	const std::string _odeApp;
	const std::string _odeTmpFile;

	void init();
	void updateInitialConditions(GrStat& stats);
	virtual void solveODE(GrStat& stats);
	void updateQueue(GrStat& stats);
	void recruitMacsGetTcellSources(GrSimulation& sim, GrStat& stats,
			ThresholdGridCellPtrList tcellSources[TCELL_TYPE_COUNT]);
	void recruitMac(GrSimulation& sim, GridCell* pSource);
	void recruitTcells(GrSimulation& sim, GrStat& stats,
			ThresholdGridCellPtrList tcellSources[TCELL_TYPE_COUNT]);

public:
	RecruitmentLnODE(const std::string& odeApp, const std::string& odeTmpFile);
	virtual ~RecruitmentLnODE();
	void recruit(GrSimulation& sim);
};

#endif /* RECRUITMENTLNODE_H_ */
