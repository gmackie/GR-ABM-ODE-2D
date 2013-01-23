#ifndef LHS_H
#define LHS_H

#include "range.h"
#include "rand.h"
extern Rand g_Rand;
#include <vector>
#include <limits>
#include <boost/property_tree/ptree.hpp>
#include <math.h>
#include <fstream>

template<typename Param, typename ParamFileHandler=typename Param::ParamFileHandler>
struct LHS {
  typedef Param param_t;
  typedef boost::property_tree::ptree ptree;
  std::vector<std::vector<int> > intBins;
  std::vector<std::vector<double> > floatBins;
  size_t nSamples;

  ptree& pt;
  Param& param;
  ParamFileHandler* handler;
  bool logScale;

  struct RangeEnumerator {
    LHS* lhs;
    RangeEnumerator(LHS* _lhs) : lhs(_lhs) {}
    template<typename T>
    inline typename boost::enable_if<boost::is_arithmetic<T> >::type
    visit(T& val, typename Param::template ParamDescriptor<T>& desc);
    template<typename T>
    inline typename boost::disable_if<boost::is_arithmetic<T> >::type
    visit(T&, typename Param::template ParamDescriptor<T>&) {}  //Skip
  };
  struct ParamBuilder {
    LHS* lhs;
    size_t intcnt, floatcnt;
    ParamBuilder(LHS* _lhs) : lhs(_lhs), intcnt(0), floatcnt(0) {}
    template<typename T>
    inline typename boost::enable_if<boost::is_integral<T> >::type
    visit(T& val, typename Param::template ParamDescriptor<T>& desc);
    template<typename T>
    inline typename boost::enable_if<boost::is_floating_point<T> >::type
    visit(T& val, typename Param::template ParamDescriptor<T>& desc);
    template<typename T>
    inline typename boost::disable_if<boost::is_arithmetic<T> >::type
    visit(T&, typename Param::template ParamDescriptor<T>&) {}  //Skip
  };

  LHS(size_t _nSamples, ptree& _pt, Param& p, ParamFileHandler* _handler, bool _logScale=false)
    : nSamples(_nSamples), pt(_pt), param(p), handler(_handler), logScale(_logScale)
  {
    assert(_nSamples > 0);
    RangeEnumerator visitor(this);
    param.visitProperties(visitor);
  }

  /// Constant int
  template<typename T>
  typename boost::enable_if<boost::is_integral<T> >::type
  addValues() { intBins.push_back(std::vector<int>()); }
  /// Constant float
  template<typename T>
  typename boost::enable_if<boost::is_floating_point<T> >::type
  addValues() { floatBins.push_back(std::vector<double>()); }
  /// Discrete int
  template<typename T>
  typename boost::enable_if<boost::is_integral<T> >::type
  addValues(std::vector<T>& values) {
    intBins.push_back(std::vector<int>(values.size()));
    std::vector<int>& bins = intBins.back();
    for(size_t i=0;i<values.size();i++) //Convert to int
      bins[i] = int(values[i]);
    for(size_t i=values.size();i<nSamples;i++)  // Add more if needed
      bins.push_back(bins[g_Rand.getInt(values.size(), 0)]);
    for(size_t i=nSamples;i<values.size();i++)  // Remove some randomly if needed
      bins.erase(bins.begin()+g_Rand.getInt(bins.size(), 0));
  }
  /// Discrete float
  template<typename T>
  typename boost::enable_if<boost::is_floating_point<T> >::type
  addValues(std::vector<T>& values) {
    floatBins.push_back(std::vector<double>(values.size()));
    std::vector<double>& bins = floatBins.back();
    for(size_t i=0;i<values.size();i++) //Convert to float
      bins[i] = double(values[i]);
    for(size_t i=values.size();i<nSamples;i++)  // Add more if needed
      bins.push_back(bins[g_Rand.getInt(values.size(), 0)]);
    for(size_t i=nSamples;i<values.size();i++)  // Remove some randomly if needed
      bins.erase(bins.begin()+g_Rand.getInt(bins.size(), 0));
  }
  /// Ranged int
  template<typename T>
  typename boost::enable_if<boost::is_integral<T> >::type
  addValues(Range<T>& range) {

    intBins.push_back(std::vector<int>(nSamples));
    std::vector<int>& bins = intBins.back();
    size_t count = nSamples / (range.max - range.min + 1);
    for(size_t j=0;j<nSamples;j++)
    {
      if(!count) {
        T sz = (range.max - range.min + 1) / nSamples;
        T min = j * sz + range.min;
        T max = std::min(range.max + 1, T(j+1) * sz + range.min);
        bins[j] = g_Rand.getInt(max, min);
      }
      else {
        T val = j / count + range.min;
        if(val <= range.max) bins[j] = int(val);
        else bins[j] = g_Rand.getInt(range.max + 1, range.min);
      }
    }
  }
  /// Ranged float
  template<typename T>
  typename boost::enable_if<boost::is_floating_point<T> >::type
  addValues(Range<T>& range) {
   floatBins.push_back(std::vector<double>(nSamples));
   std::vector<double>& bins = floatBins.back();
   if(logScale) 
     for(size_t i=0;i<nSamples;i++)
     {
        const double min = log10(range.min);
        const double max = log10(range.max);
        const double a = min + i * (max - min) / nSamples;
        const double b = min + (i+1) * (max - min) / nSamples;
        bins[i] = pow(double(10), double(g_Rand.getReal(a, b)));
     }
   else
     for(size_t i=0;i<nSamples;i++)
     {
        const double a = range.min + i * (range.max - range.min) / nSamples;
        const double b = range.min + (i+1) * (range.max - range.min) / nSamples;
        bins[i] = double(g_Rand.getReal(a, b));
     }
  }
  void performLHS();
};

template<typename Param, typename ParamFileHandler>
template<typename T>
inline typename boost::enable_if<boost::is_arithmetic<T> >::type
LHS<Param, ParamFileHandler>::RangeEnumerator::visit(T& val, typename Param::template ParamDescriptor<T>& desc) {
  //Get the range value from pt
  boost::optional<ptree&> kid = lhs->pt.get_child_optional(lhs->handler->getFullPath(desc.name, desc.path));
  if(!!kid && kid->data()[0] == '[') { //It's a range
    Range<T> range = kid->get_value<Range<T> >();
    if(!!desc.range && !desc.range->contains(range))
      throw std::runtime_error(std::string("Invalid range for parameter: ") + desc.path+'.'+desc.name);

    desc.wasRead = true;

    if(range.min == range.max)
    {
      val = range.min;
    }
    else { //Build the sample list
      lhs->addValues(range);
      return;
    }
  }
  else if(!!kid) {
    std::stringstream ss(kid->data());
    std::istream_iterator<T> begin(ss), end;
    std::vector<T> values(begin, end);
    if(values.size() == 0)
      throw std::runtime_error(std::string("Invalid format for parameter: ")+desc.path+'.'+desc.name);
    else if(values.size() == 1)
    {
      val = values[0];  //Sample list is constant
      desc.wasRead = true;
    }
    else {
      if(!!desc.range)
        for(typeof(values.begin()) it = values.begin();it!=values.end();it++)
          if(!desc.range->contains(*it))
            throw std::runtime_error(std::string("Invalid value for ") + desc.path+'.'+desc.name+": "+ boost::lexical_cast<std::string>(*it));
      //Build the sample list
      lhs->addValues(values);
      desc.wasRead = true;
      return;
    }
  }
  else if(!!desc.def)
    {
      val = *(desc.def);
    }
  else {
    std::cerr<<"Warning - Parameter value cannot be determined: "<<desc.path<<'.'<<desc.name<<std::endl;
    val = 0;
  }
  lhs->addValues<T>();  //Determined to be constant
}

/// Picks a bin and removes it
template<typename Param, typename ParamFileHandler>
template<typename T>
inline typename boost::enable_if<boost::is_integral<T> >::type
LHS<Param, ParamFileHandler>::ParamBuilder::visit(T& val, typename Param::template ParamDescriptor<T>&) {
  assert(intcnt < lhs->intBins.size());
  std::vector<int>& bins = lhs->intBins[intcnt++];
  if(bins.size() == 0) return;  // Constant value
  std::vector<int>::iterator it = bins.begin() + g_Rand.getInt(bins.size(), 0);
  val = T(*it);
  bins.erase(it);
}

/// Picks a bin and removes it
template<typename Param, typename ParamFileHandler>
template<typename T>
inline typename boost::enable_if<boost::is_floating_point<T> >::type
LHS<Param, ParamFileHandler>::ParamBuilder::visit(T& val, typename Param::template ParamDescriptor<T>&) {
  assert(floatcnt < lhs->floatBins.size());
  std::vector<double>& bins = lhs->floatBins[floatcnt++];
  if(bins.size() == 0) return;  // Constant value
  std::vector<double>::iterator it = bins.begin() + g_Rand.getInt(bins.size(), 0);
  val = T(*it);
  bins.erase(it);
}

template<typename Param, typename ParamFileHandler>
inline void LHS<Param, ParamFileHandler>::performLHS() {
  for(size_t i=0;i<nSamples;i++) 
  {
    ParamBuilder visitor(this);
    param.visitProperties(visitor);

    if (!param.validate())
    {
      std::cerr << "Sample " << i << " had validation errors." << std::endl;
      exit(1);
    }

    std::ofstream file((boost::lexical_cast<std::string>(i+1) + ".xml").c_str());
    param.save(file, handler, pt);
  }
}

#endif // LHS_H
