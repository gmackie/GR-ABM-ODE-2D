/*
 * areatest.cpp
 *
 *  Created on: Jan 27, 2010
 *      Author: mohammed
 */

#include "areatest.h"

AreaTest::AreaTest(double alpha, int testPeriod, int samplePeriod)
	: TTest(alpha, testPeriod, samplePeriod, "OUTCOME_AREA")
{
}

AreaTest::~AreaTest()
{
}

void AreaTest::update(const int time, const int index, Stats& stats)
{
	TTest::update(time, index, stats, stats.getAreaTNF());
}

void AreaTest::evaluate(const int index, Stats& stats, double degressOfFreedom)
{

	//double pValContainment = twoSided(degressOfFreedom);
	//double pValClearance = oneSidedMore(degressOfFreedom);
	double pValDissemination = oneSidedLess(degressOfFreedom);

	//printf("Containment: %lf\tClearance: %lf\tDissemination: %lf\n",
	//		pValContainment, pValClearance, pValDissemination);

	/*else if (pValContainment > _alpha)
	{
		stats.setGrStatus(GR_CONTAINMENT);
	}
	else if (pValClearance < _alpha)
	{
		stats.setGrStatus(GR_CLEARANCE);
	}*/
	if (pValDissemination < _alpha)
	{
		stats.setStatus(index, GR_DISSEMINATION);
	}
	else //if (_stat1.getMean() == _stat2.getMean() && _stat1.getVariance() == _stat2.getVariance())
	{
		stats.setStatus(index, GR_CONTAINMENT);
	}
	/*else
	{
		stats.setGrStatus(GR_UNKNOWN);
	}*/

}
