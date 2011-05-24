/*
 * rand.h
 *
 *  Created on: 06-nov-2009
 *      Author: M. El-Kebir
 */

#ifndef RAND_H
#define RAND_H

#include <assert.h>
#include <iostream>

#ifdef __APPLE__
	#include <boost/random.hpp>

	typedef boost::mt19937 mt19937_t;
	typedef boost::random::ranlux64_base_01 ranlux64_base_01_t;
	typedef boost::uniform_real<double> uniform_real_t;
	typedef boost::uniform_int<int> uniform_int_t;
#else
	#ifdef _MSC_VER
		#include <random>
	#else
		#include <boost/tr1/random.hpp>
	#endif

	typedef std::tr1::mt19937 mt19937_t;
	typedef std::tr1::ranlux64_base_01 ranlux64_base_01_t;
	typedef std::tr1::uniform_real<double> uniform_real_t;
	typedef std::tr1::uniform_int<int> uniform_int_t;
#endif

class Rand
{
private:
	/*
	 * !!! If the data members change then the serialize and deserialize functions need to be updated !!!
	 */
	unsigned int _seed;
	mt19937_t _eng;
	ranlux64_base_01_t _realEng;
	uniform_real_t _uniformDist01;


public:
	Rand(unsigned int seed);
	~Rand();
	void setSeed(unsigned int seed);
	unsigned int getSeed();

	// Needed so we can use this class with the random_shuffle function in the C++ algorithms library.
	ptrdiff_t operator() (ptrdiff_t max);

	double getReal();
	double getReal(double a, double b);
	int getInt(int b, int a = 0);
	void serialize(std::ostream& out) const;
	void deserialize(std::istream& in);
	std::size_t getSerialSize() const;
	void test(int time);
};

inline void Rand::setSeed(unsigned int seed)
{
	_seed = seed;
	_eng.seed(seed);
	_realEng.seed(seed);
}

inline unsigned int Rand::getSeed()
{
	return _seed;
}

inline ptrdiff_t Rand::operator() (ptrdiff_t max)
{
	return getInt(max, 0);
}

inline double Rand::getReal()
{
	return _uniformDist01(_realEng);
}

inline double Rand::getReal(double a, double b)
{
	uniform_real_t rnd(a, b);
	return rnd(_realEng);
}

inline int Rand::getInt(int b, int a)
{
	assert(a < b);

	uniform_int_t rnd(a, b - 1);
	return rnd(_eng);
}

inline std::size_t Rand::getSerialSize() const
{
	// We don't serialize _uniformDist01
	return sizeof(Rand) - sizeof(uniform_real_t);
}

#endif /* RAND_H */
