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
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/split_free.hpp>

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
  static const std::string _ClassName;

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
  double getLogNormal(double mean, double sigma);
  void test(int time);
  template<class Archive>
  void serialize(Archive& ar, const unsigned int version);
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

inline double Rand::getLogNormal(double mean, double sigma)
{
  boost::lognormal_distribution<double> rnd(mean, sigma);
  return rnd(_realEng);
}

namespace boost { namespace serialization {
  /*
   * Some serialization methods defined because boost didn't integrate random
   * library with serialization library.  Basically splits the serialization
   * process in order to use the iostream version which *is* defined
   */
  template<class Archive>
  inline void serialize(Archive& ar, ranlux64_base_01_t& realEng, const unsigned int version)
  {
    split_free(ar, realEng, version);
  }
  template<class Archive>
  inline void save(Archive& ar, const ranlux64_base_01_t& realEng, const unsigned int /*version*/)
  {
    std::ostringstream ss;
    ss<<realEng;
    std::string s(ss.str());
    ar << boost::serialization::make_nvp("realEng", s);
  }
  template<class Archive>
  inline void load(Archive& ar, ranlux64_base_01_t& realEng, const unsigned int /*version*/)
  {
    std::string s;
    ar >> boost::serialization::make_nvp("realEng", s);
    std::istringstream ss(s);
    ss >> realEng;
  }
  template<class Archive>
  inline void serialize(Archive& ar, mt19937_t& eng, const unsigned int version)
  {
    split_free(ar, eng, version);
  }
  template<class Archive>
  inline void save(Archive& ar, const mt19937_t& eng, const unsigned int /*version*/)
  {
    std::ostringstream ss;
    ss << eng;
    std::string s(ss.str());
    ar << boost::serialization::make_nvp("eng", s);
  }
  template<class Archive>
  inline void load(Archive& ar, mt19937_t& eng, const unsigned int /*version*/)
  {
    std::string s;
    ar >> boost::serialization::make_nvp("eng", s);
    std::istringstream ss(s);
    ss >> eng;
  }
}}

template<class Archive>
void Rand::serialize(Archive& ar, const unsigned int /*version*/)
{
  ar & BOOST_SERIALIZATION_NVP(_seed);
  ar & BOOST_SERIALIZATION_NVP(_eng);
  ar & BOOST_SERIALIZATION_NVP(_realEng);
}

#endif /* RAND_H */
