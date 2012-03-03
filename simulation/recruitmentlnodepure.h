/*
 * RecruitmentLnODEPure.h
 *
 *  Created on: Jul 4, 2010
 *      Author: mohammed
 */

#ifndef RECRUITMENTLNODEPURE_H_
#define RECRUITMENTLNODEPURE_H_

#include "recruitmentlnode.h"

class RecruitmentLnODEPure : public RecruitmentLnODE
{
private:
	static const int _idxTH0 = 0;
	static const int _idxTH1 = 1;
	static const int _idxT80 = 2;
	static const int _idxT8 = 3;
	static const int _idxPrecursorCTL = 4;
	static const int _idxEffectorTH0 = 5;
	static const int _idxEffectorT80 = 6;

	virtual void solveODE(const int time, GrStat& stats);

public:
	RecruitmentLnODEPure();
	virtual ~RecruitmentLnODEPure();

};

#endif /* RECRUITMENTLNODEPURE_H_ */
