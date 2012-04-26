/*
 * ttest.cpp
 *
 *  Created on: Feb 4, 2010
 *      Author: mohammed
 */

#include "ttest.h"
#include <boost/math/distributions/students_t.hpp>

TTest::TTest(double alpha, int testPeriod, int samplePeriod, std::string methodName)
	: _stat1()
	, _stat2()
	, _alpha(alpha)
	, _testPeriod(testPeriod)
	, _samplePeriod(samplePeriod)
	, _methodName(methodName)
{
}

TTest::~TTest()
{
}

double TTest::degreesOfFreedom() const
{
	double v = _stat1.getVariance() / _stat1.getNrSamples() +
		_stat2.getVariance() / _stat2.getNrSamples();

	v *= v;

	double t1 = _stat1.getVariance() / _stat1.getNrSamples() + DBL_EPSILON;

	t1 *= t1;
	t1 /= (_stat1.getNrSamples() + 1);

	double t2 = _stat2.getVariance() / _stat2.getNrSamples() + DBL_EPSILON;

	t2 *= t2;
	t2 /= (_stat2.getNrSamples() + 1);

	v = v / (t1 + t2) - 2;

	//printf("v = %lf\n", v);
	//printf("n_1 = %d, n_2 = %d\n", _stat1.getNrSamples(), _stat2.getNrSamples());
	//printf("\\mu_1 = %lf, \\mu_2 = %lf\n", _stat1.getMean(), _stat2.getMean());
	//printf("\\sigma^2_1 = %lf, \\sigma^2_2 = %lf\n", _stat1.getVariance(), _stat2.getVariance());

	return v < 0 ? 10 * DBL_EPSILON : v;
}

double TTest::tStat() const
{
	double t = _stat1.getMean() - _stat2.getMean();

	t /= sqrt(_stat1.getVariance() / _stat1.getNrSamples() +
			_stat2.getVariance() / _stat2.getNrSamples());

	return t;
}

double TTest::twoSided(double v) const
{
	using namespace boost::math;
	students_t dist(v);

	return 2 * cdf(complement(dist, fabs(tStat())));
}

double TTest::oneSidedLess(double v) const
{
	using namespace boost::math;
	students_t dist(v);

	return cdf(dist, tStat());
}

double TTest::oneSidedMore(double v) const
{
	using namespace boost::math;
	students_t dist(v);

	return cdf(complement(dist, tStat()));
}

void TTest::perform(const int index, Stats& stats)
{

	// This can happen when we do one last t-test at the end of
	// a simulation and we've already done one recently enough
	// that we don't have any samples for this one.
	// If there is only one sample should we skip that also?
	if (_stat2.getNrSamples() == 0)
	{
		return;
	}

	if (_stat1.getNrSamples() == 0)
	{
		_stat1 = _stat2;
	}

	// If no bacteria it always clearance, regardless of the type of t-test used.
	if (stats.getTotExtMtb() + stats.getTotIntMtb() == 0)
	{
		stats.setStatus(index, GR_CLEARANCE);
		return;
	}

	double v = degreesOfFreedom();

	evaluate(index, stats, v);

	_stat1 = _stat2;
	_stat2.reset();
}

void TTest::update(const int time, const int index, Stats& stats, double value)
{
	if (time % _samplePeriod == 0)
	{
		_stat2.update(value);
	}

	if (time % _testPeriod == 0 && _stat2.getNrSamples() > 1)
	{
		perform(index, stats);
	}
}

