/*
 * mtbtest.h
 *
 *  Created on: Feb 5, 2010
 *      Author: mohammed
 */

#ifndef MTBTEST_H_
#define MTBTEST_H_

#include "ttest.h"

class MtbTest: public TTest
{
public:
	MtbTest(double alpha, int testPeriod, int samplePeriod);
	~MtbTest();
	void update(const int time, const int index, GrStat& stats);
	void evaluate(const int index, GrStat& stats, double degressOfFreedom);
	OutcomeMethod getMethod() const;

};

inline OutcomeMethod MtbTest::getMethod() const
{
	return OUTCOME_MTB;
}

#endif /* MTBTEST_H_ */
