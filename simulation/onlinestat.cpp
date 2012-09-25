/*
 * onlinestat.cpp
 *
 *  Created on: Jan 27, 2010
 *      Author: mohammed
 */

#include "onlinestat.h"

OnlineStat::OnlineStat()
  : _n(0)
  , _mean(0)
  , _m2(0)
{
}

OnlineStat::~OnlineStat()
{
}

void OnlineStat::reset()
{
  _n = 0;
  _mean = 0;
  _m2 = 0;
}

void OnlineStat::update(double x)
{
  double delta = x - _mean;

  _n++;
  _mean = _mean + delta / _n;
  _m2 = _m2 + delta * (x - _mean);
}
