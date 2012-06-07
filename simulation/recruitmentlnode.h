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
#include "stat.h"
#include "recruitmentbase.h"
#include <string>

class RecruitmentLnODE : public RecruitmentBase
{
private:
	static const std::string _ClassName;

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
		TcellTypePair(int birthtime, TcellType type) : _birthtime(birthtime), _type(type) {};
		int _birthtime;
		TcellType _type;
	};

	std::vector<TcellTypePair> _tcellQueue;

	int _prevMiMci;
	double _odeInitialConditions[_nrConditions];
	const std::string _odeApp;
	const std::string _odeTmpFile;

	void init();
	void updateInitialConditions(Stats& stats);
	virtual void solveODE(const int time, const Stats& statsPrevious, Stats& stats);
	void updateQueue(const int time, Stats& stats);
	void recruitMacsGetTcellSources(GrSimulation& sim, Stats& stats,
			ThresholdPosList tcellSources[TCELL_TYPE_COUNT]);
	void recruitMac(GrSimulation& sim, const Pos& pSource);
	void recruitTcells(GrSimulation& sim, Stats& stats,
			ThresholdPosList tcellSources[TCELL_TYPE_COUNT]);

public:
	RecruitmentLnODE(const std::string& odeApp, const std::string& odeTmpFile, std::istream& in);
	RecruitmentLnODE(const std::string& odeApp, const std::string& odeTmpFile);
	virtual ~RecruitmentLnODE();

	RecruitmentMethod getMethod() const;
	void serialize(std::ostream& out) const;
	void deserialize(std::istream& in);

	void recruit(GrSimulation& sim);
};

inline RecruitmentMethod RecruitmentLnODE::getMethod() const
{
	return RECR_LN_ODE;
}

#endif /* RECRUITMENTLNODE_H_ */
