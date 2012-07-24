/*
 * vector.h
 *
 *  Created on: 21-sep-2008
 *      Author: S030858
 */

#ifndef VECTOR_H_
#define VECTOR_H_

#include <assert.h>
#include <math.h>

namespace math
{

template<int N = 3, typename T = float>
class vector
{
private:
  T _v[N];

public:
  vector();
  vector(T x);
  vector(const vector<N, T>& vec);
  vector(T(&a)[N]);
  ~vector();
  const T* v() const;
  T operator [](int i) const;
  T& operator [](int i);
  void mult(T c);
  vector& operator *=(T c);
  void add(const vector<N, T>& vec);
  vector& operator +=(const vector<N, T>& vec);
  void min(const vector<N, T>& vec);
  vector& operator -=(const vector<N, T>& vec);
  float magnitude() const;
  void normalize();
  double dotproduct(const vector<N, T>& vec) const;
};

template<int N, typename T>
inline vector<N, T>::vector()
{
  for (int i = 0; i < N; i++)
    {
      _v[i] = 0;
    }
}

template<int N, typename T>
inline vector<N, T>::vector(T x)
{
  for (int i = 0; i < N; i++)
    {
      _v[i] = x;
    }
}

template<int N, typename T>
inline vector<N, T>::vector(const vector<N, T>& vec)
{
  for (int i = 0; i < N; i++)
    {
      this->_v[i] = vec._v[i];
    }
}

template<int N, typename T>
inline vector<N, T>::vector(T(&a)[N])
{
  for (int i = 0; i < N; i++)
    {
      _v[i] = a[i];
    }
}

template<int N, typename T>
inline vector<N, T>::~vector()
{
}

template<int N, typename T>
inline const T* vector<N, T>::v() const
{
  return _v;
}

template<int N, typename T>
inline T vector<N, T>::operator [](int i) const
{
  assert(0 <= i && i < N);
  return _v[i];
}

template<int N, typename T>
inline T& vector<N,T>::operator [](int i)
{
  assert(0 <= i && i < N);
  return _v[i];
}

template<int N, typename T>
inline void vector<N, T>::mult(T c)
{
  for (int i = 0; i < N; i++)
    {
      _v[i] *= c;
    }
}

template<int N, typename T>
inline vector<N, T>& vector<N, T>::operator *=(T c)
{
  this->mult(c);
  return *this;
}

template<int N, typename T>
inline void vector<N, T>::add(const vector<N, T>& vec)
{
  for (int i = 0; i < N; i++)
    {
      _v[i] += vec[i];
    }
}

template<int N, typename T>
inline vector<N, T>& vector<N, T>::operator +=(const vector<N, T>& vec)
{
  this->add(vec);
  return *this;
}

template<int N, typename T>
inline void vector<N, T>::min(const vector<N, T>& vec)
{
  for (int i = 0; i < N; i++)
    {
      _v[i] -= vec[i];
    }
}

template<int N, typename T>
inline vector<N, T>& vector<N, T>::operator -=(const vector<N, T>& vec)
{
  this->min(vec);
  return *this;
}

template<int N, typename T>
inline float vector<N, T>::magnitude() const
{
  float sum = 0.0;
  for (int i = 0; i < N; i++)
    {
      sum += (_v[i] * _v[i]);
    }

  return sqrt(sum);
}

template<int N, typename T>
inline double vector<N, T>::dotproduct(const vector<N, T>& vec) const
{
  double sum = 0;
  for (int i = 0; i < N; i++)
    {
      sum += ((double) vec._v[i] * (double) _v[i]);
    }
  return sum;
}

template<int N, typename T>
inline void vector<N, T>::normalize()
{
  float mag = magnitude();
  for (int i = 0; i < N; i++)
    {
      _v[i] /= mag;
    }
}

template<int N, typename T>
vector<N, T> operator*(T c, const vector<N, T>& vec)
{
  vector<N, T> res = vec;
  return res *= c;
}

template<int N, typename T>
vector<N, T> operator*(const vector<N, T>& vec, T c)
{
  vector<N, T> res = vec;
  return res *= c;
}

template<int N, typename T>
vector<N, T> operator+(const vector<N, T>& vec1, const vector<N, T>& vec2)
{
  vector<N, T> res = vec1;
  return res += vec2;
}

template<int N, typename T>
vector<N, T> operator-(const vector<N, T>& vec1, const vector<N, T>& vec2)
{
  vector<N, T> res = vec1;
  return res -= vec2;
}

}

typedef math::vector<2, float> vec2f;
typedef math::vector<3, float> vec3f;
typedef math::vector<2, int> vec2i;
typedef math::vector<3, int> vec3i;

#endif /* VECTOR_H_ */
