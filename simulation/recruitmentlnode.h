/*
 * recruitmentlnode.h
 *
 *  Created on: Jun 17, 2010
 *      Author: mohammed
 */

/*
 * This class originally was used to run a lymph node ode that used a DLL created by Matlab.
 * The Matlab DLL is no longer used, but everything in this function is used by subclasses,
 * except that they override the solveODE function.
 *
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

	static const int _nrConditions = 43;

	static const int _idxMDC = 13;
	static const int _idxNaiveCD4 = 14;
	static const int _idxNaiveCD8 = 16;
	static const int _idxEffectorTH1 = 39;
	static const int _idxEffectorT8 = 41;
	static const int _idxCTL = 42;

	double _tcellTable[TCELL_TYPE_COUNT]; // contains lower bounds

	int _tcellQueueCount[TCELL_TYPE_COUNT];
	struct TcellTypePair
	{
		int _birthtime;
		TcellType _type;
	};

	std::vector<TcellTypePair> _tcellQueue;

	int _prevMiMci;
	double _odeInitialConditions[_nrConditions];
	const std::string _odeApp;
	const std::string _odeTmpFile;

	void init();
	void updateInitialConditions(GrStat& stats);
	virtual void solveODE(const int time, const GrStat& statsPrevious, GrStat& stats);
	void updateQueue(const int time, GrStat& stats);
	void recruitMacsGetTcellSources(GrSimulation& sim, GrStat& stats,
			ThresholdPosList tcellSources[TCELL_TYPE_COUNT]);
	void recruitMac(GrSimulation& sim, const Pos& pSource);
	void recruitTcells(GrSimulation& sim, GrStat& stats,
			ThresholdPosList tcellSources[TCELL_TYPE_COUNT]);

public:
	RecruitmentLnODE(const std::string& odeApp, const std::string& odeTmpFile);
	virtual ~RecruitmentLnODE();
	void recruit(GrSimulation& sim);
};

#endif /* RECRUITMENTLNODE_H_ */
