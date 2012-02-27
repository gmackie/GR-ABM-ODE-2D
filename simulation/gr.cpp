/*
 * gr.cpp
 *
 *  Created on: 03-nov-2009
 *      Author: M. El-Kebir
 */

#include "gr.h"
#include <sys/time.h>

Rand g_Rand = Rand(1);

// Create a seed based on the current time.
// If we just use time in seconds we will get duplicate seeds when doing a large number of runs,
// such as running a big LHS on a cluster.
// If we use just the micro-second part of the high resolution time we only get seeds that vary
// between 1 and 1,000,000.
// So we XOR the two times, which gives a larger range of seeds, with the lower order part varying
// with a higher resolution than we would get just by using seconds.
unsigned int createTimeSeed()
{
	 struct timeval curTimeHiRes;
	 if (gettimeofday(&curTimeHiRes, NULL) != 0)
			throw std::runtime_error("Error getting current high resolution time.");

	 unsigned int seed = curTimeHiRes.tv_sec ^ curTimeHiRes.tv_usec;

	return seed;
}
