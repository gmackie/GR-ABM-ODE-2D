/*
 * areatest.h
 *
 *  Created on: Jan 27, 2010
 *      Author: mohammed
 */

#ifndef AREATEST_H_
#define AREATEST_H_

#include "onlinestat.h"
#include "stat.h"
#include "ttest.h"

class AreaTest : public TTest
{
public:
	AreaTest(double alpha, int testPeriod, int samplePeriod);
	~AreaTest();
	void update(const int time, const int index, Stats& stats);
	void evaluate(const int index, Stats& stats, double degressOfFreedom);
	OutcomeMethod getMethod() const;
};

inline OutcomeMethod AreaTest::getMethod() const
{
	return OUTCOME_AREA;
}

#endif /* AREATEST_H_ */
