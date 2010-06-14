/*
 * onlinestat.h
 *
 *  Created on: Jan 27, 2010
 *      Author: mohammed
 */

#ifndef ONLINESTAT_H_
#define ONLINESTAT_H_

#include "gr.h"

/* Algorithm is due to Donald E. Knuth (1998).
 * The Art of Computer Programming,
 * volume 2: Seminumerical Algorithms,
 * 3rd edition, p. 232. Boston: Addison-Wesley
 */
class OnlineStat
{
private:
	int _n;
	double _mean;
	double _m2;

public:
	OnlineStat();
	virtual ~OnlineStat();
	void reset();
	void update(double x);
	int getNrSamples() const;
	double getMean() const;
	double getVariance() const;
};

inline int OnlineStat::getNrSamples() const
{
	return _n;
}

inline double OnlineStat::getMean() const
{
	return _mean;
}

inline double OnlineStat::getVariance() const
{
	return _m2 / (_n - 1);
}

#endif /* ONLINESTAT_H_ */
