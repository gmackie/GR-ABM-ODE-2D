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
  void update(const int time, const int index, Stats& stats);
  void evaluate(const int index, Stats& stats, double degressOfFreedom);
  OutcomeMethod getMethod() const;
  /*virtual*/
  TTest* clone() const
  {
    return new MtbTest(*this);
  }
};

inline OutcomeMethod MtbTest::getMethod() const
{
  return OUTCOME_MTB;
}

#endif /* MTBTEST_H_ */
