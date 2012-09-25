/*
 * ttest.h
 *
 *  Created on: Feb 4, 2010
 *      Author: mohammed
 */

#ifndef TTEST_H_
#define TTEST_H_

#include "onlinestat.h"
#include "stat.h"

/* This class provides functionality
 * for performing Welch's Student-t test */
class TTest
{
protected:
  OnlineStat _stat1;
  OnlineStat _stat2;
  double _alpha;
  int _testPeriod;
  int _samplePeriod;
  std::string _methodName;

  void update(const int time, const int index, Stats& stats, double value);
  double degreesOfFreedom() const;
  double tStat() const;
  // returns p-value; H_1 : \mu_1 != \mu_2
  double twoSided(double v) const;
  // returns p-value; H_1 : \mu_1 < \mu_2
  double oneSidedLess(double v) const;
  // returns p-value; H_1 : \mu_1 > \mu_2
  double oneSidedMore(double v) const;
  virtual void evaluate(const int index, Stats& stats, double degressOfFreedom) = 0;

public:
  TTest(double alpha, int testPeriod, int samplePeriod, std::string methodName);
  virtual ~TTest();
  virtual void update(const int time, const int index, Stats& stats) = 0;
  void perform(const int index, Stats& stats);
  virtual OutcomeMethod getMethod() const = 0;
  virtual TTest* clone() const = 0;
  std::string getMethodName() const;
  void setAlpha(double alpha);
  void setTestPeriod(int testPeriod);
  void setSamplePeriod(int samplePeriod);
  double getAlpha() const;
  int getTestPeriod() const;
  int getSamplePeriod() const;
};

inline double TTest::getAlpha() const
{
  return _alpha;
}

inline int TTest::getTestPeriod() const
{
  return _testPeriod;
}

inline int TTest::getSamplePeriod() const
{
  return _samplePeriod;
}

inline std::string TTest::getMethodName() const
{
  return _methodName;
}


inline void TTest::setAlpha(double alpha)
{
  _alpha = alpha;
}

inline void TTest::setTestPeriod(int testPeriod)
{
  _testPeriod = testPeriod;
}

inline void TTest::setSamplePeriod(int samplePeriod)
{
  _samplePeriod = samplePeriod;
}

#endif /* TTEST_H_ */
