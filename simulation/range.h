#ifndef RANGE_H
#define RANGE_H

#include <iostream>
#include <assert.h>

/// Class to hold valid ranges for parameters
template<typename T>
struct Range {
  typedef T value_type;
  T min, max;
  Range(const T& _min = T()) : min(_min), max(_min) {}
  Range(const T& _min, const T& _max) : min(std::min(_min, _max)), max(std::max(_min, _max)) {}
  template<typename U> Range(const Range<U>& o) : min(o.min), max(o.max) {}
  template<typename U> operator Range<U>() { return Range<U>(min, max);  }  //Should compile error if T isn't convertible to U
  bool contains(const T& v) const { return min<=v && v<=max; }
  template<typename U> bool contains(const Range<U>& v) const { return contains(v.min) && contains(v.max); }
};
template<typename T>
inline std::ostream& operator<<(std::ostream& s, Range<T>& r) {
  return s<<'['<<r.min<<','<<r.max<<']';
}

template<typename T>
inline std::istream& operator>>(std::istream& s, Range<T>& r) {
  char c;
  s>>c; assert(c=='[');
  s>>r.min;
  s>>c; assert(c==',');
  s>>r.max;
  s>>c; assert(c==']');
  return s;
}

#endif // RANGE_H
