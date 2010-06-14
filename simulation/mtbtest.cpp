/*
 * mtbtest.cpp
 *
 *  Created on: Feb 5, 2010
 *      Author: mohammed
 */

#include "mtbtest.h"

MtbTest::MtbTest(double alpha, int testPeriod, int samplePeriod)
	: TTest(alpha, testPeriod, samplePeriod, "OUTCOME_MTB")
{
}

MtbTest::~MtbTest()
{
}

void MtbTest::update(const int time, const int index, GrStat& stats)
{
	TTest::update(time, index, stats, stats.getTotExtMtb() + stats.getTotIntMtb());
}

void MtbTest::evaluate(const int index, GrStat& stats, double degressOfFreedom)
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
		if (stats.getTotExtMtb() >= stats.getTotIntMtb())
		{
			stats.setGrStatus(index, GR_DISSEMINATION);
		}
		else
		{
			stats.setGrStatus(index, GR_DISSEMINATION_INCONSISTENT);
		}
	}
	else //if (_stat1.getMean() == _stat2.getMean() && _stat1.getVariance() == _stat2.getVariance())
	{
		if (stats.getTotExtMtb() <= stats.getTotIntMtb())
		{
			stats.setGrStatus(index, GR_CONTAINMENT);
		}
		else
		{
			stats.setGrStatus(index, GR_CONTAINMENT_INCONSISTENT);
		}
	}
	/*else
	{
		stats.setGrStatus(GR_UNKNOWN);
	}*/


}


