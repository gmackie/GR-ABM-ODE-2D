#ifndef STAT_H
#define STAT_H
//Heavy header file, include only in cpp files, forward declare otherwise
#include <boost/mpl/find.hpp>
#include <boost/array.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/preprocessor/expr_if.hpp>
#include <boost/preprocessor/comparison/equal.hpp>
#include <iostream>
#include <numeric>
#include "gr.h"
#include "agent.h"
#include "macrophage.h"
#include "tregulatory.h"
#include "tcytotoxic.h"
#include "tgamma.h"

typedef enum {	GR_CONTAINMENT, GR_CONTAINMENT_INCONSISTENT, GR_CLEARANCE,
				GR_DISSEMINATION, GR_DISSEMINATION_INCONSISTENT, GR_UNKNOWN, GR_NONE} GrStatus;
namespace ba = boost::accumulators;


// --- iostream members for boost::array<> ---
inline std::ostream& operator<<(std::ostream& s, const GrStatus& a) {
  s<<(int)a;
  return s;
}
inline std::istream& operator>>(std::istream& s, GrStatus& a) {
  int x;
  s>>x;
  a = (GrStatus)x;
  return s;
}
template<typename T, size_t N>
inline std::ostream& operator<<(std::ostream& s, const boost::array<T, N>& a) {
  s<<N<<' ';
  for(typename boost::array<T, N>::const_iterator it = a.begin(); it != a.end(); it++)
    s<<*it<<' ';
  return s;
}

template<typename T, size_t N>
inline std::istream& operator>>(std::istream& s, boost::array<T, N>& a) {
  size_t n=0;
  s>>n;
  if(n!=N)
    throw std::ios_base::failure("Array dimension misread");
  for(typename boost::array<T, N>::iterator it = a.begin(); it != a.end(); it++)
  {
    T& v = *it;
    s>>v;
  }
  return s;
}

#define RESET_AGENT_STATS
#define RESET_REGULAR_STATS

class Stats
{
  friend class boost::serialization::access;
public:
  typedef ba::accumulator_set<Scalar, ba::stats< ba::features< ba::tag::variance, ba::tag::min, ba::tag::max, ba::tag::median > > > Stat;
protected:
// --- Data Members ---
#define STAT(type, name, desc, reset) \
  type _##name;
#define ARRAY_STAT(type, name, sz, desc, reset) \
  boost::array<type, sz> _##name;
#define AGENT_STAT(type, name, desc, reset) \
  ARRAY_STAT(type, name, NAGENTS, desc, reset)
#define STATE_STAT(type, name, agent_t, desc, reset) \
  ARRAY_STAT(type, name, agent_t::NSTATES, desc, reset)
/*#define AGENT_STAT_MAKER(r, name, agent_t) TYPE _ ## agent_t ## name [agent_t::NSTATES];
///TODO: Implement AGENTSTATE_STAT
#define AGENTSTATE_STAT(type, name, desc, reset) \
  AGENT_STAT(type, Agent ## name, desc, reset)  \
  #define _TYPE type \
  BOOST_PP_SEQ_FOR_EACH(AGENT_STAT_MAKER, name, AGENT_SEQ)  \
  #undef _TYPE */
#include "stat.def"

//--- Custom Members ---
  boost::array<Stat, Mac::NSTATES> _macIntMtbStats;
	std::vector<unsigned> _intMtbFreq;

public:

// --- Constructor ---

Stats() :
	_intMtbFreq(int(_PARAM(PARAM_MAC_THRESHOLD_BURST_CI_INTMTB)) + 1, 0)
{
  clear();
}

// --- Accessors ---

#define STAT(type, name, desc, reset) \
  const type& get ## name () const { return _ ## name; }  \
  type& get ## name () { return _ ## name; }  \
  void set ## name (const type& v) { _ ## name = v; }

#define ARRAY_STAT(type, name, sz, desc, reset) \
  const type& get ## name (size_t i) const { assert(i < sz); return _ ## name [i]; }  \
  type& get ## name (size_t i) { assert(i < sz); return _ ## name [i]; }  \
  const boost::array<type, sz>& get ## name ## Array () const { return _ ## name; }  \
  void set ## name (size_t i, const type& v) { assert(i < sz); (_ ## name [i]) = v; } \
  void set ## name (const type& v) { std::fill_n((_ ## name).begin(), (size_t)sz, v); }

#define STATE_STAT(type, name, agent_t, desc, reset) \
  const type& get ## name (agent_t::State i) const { return _ ## name [i]; }  \
  type& get ## name (agent_t::State i) { return _ ## name [i]; }  \
  type get ## name () const { return std::accumulate(_ ## name .begin(), _ ## name .end(), type()); } \
  const boost::array<type, agent_t::NSTATES>& get ## name ## Array () const { return _ ## name; }  \
  void set ## name (agent_t::State i, const type& v) { (_ ## name [i]) = v; } \
  void set ## name (const type& v) { std::fill_n((_ ## name).begin(), (size_t)agent_t::NSTATES, v); }
/*  ///TODO: Implement AGENTSTATE_STAT
#define AGENTSTATE_STAT(type, name, sz, desc, reset) \
  template<typename agent_t>  \
  void inc ## name (typename agent_t::State i) { ++ _ ## name [i]; }  \
  template<typename agent_t>  \
  void dec ## name (typename agent_t::State i) { -- _ ## name [i]; }  \
  template<typename agent_t>  \
  const type& get ## name (typename agent_t::State i) const { return _ ## name [i]; }  \
  template<typename agent_t>  \
  type& get ## name (typename agent_t::State i) { return _ ## name [i]; }  \
  template<typename agent_t>  \
  void set ## name (typename agent_t::State i, const type& v) { _ ## name [i] = v; }
*/
#define AGENT_STAT(t, name, desc, reset) \
  template<typename agent_t>  \
  const t& get ## name () const { return _ ## name [boost::mpl::find<AgentTypes, agent_t>::type::pos::value]; }  \
  template<typename agent_t>  \
  t& get ## name () { return _ ## name [boost::mpl::find<AgentTypes, agent_t>::type::pos::value]; }  \
  const boost::array<t, NAGENTS>& get ## name ## Array () const { return _ ## name; }  \
  t getTot ## name () const { return std::accumulate(_ ## name .begin(), _ ## name .end(), t()); } \
  template<typename agent_t>  \
  void set ## name (const t& v) { _ ## name [boost::mpl::find<AgentTypes, agent_t>::type::pos::value] = v; } \
  /* Non-Template enum versions */ \
  const t& get ## name (AgentType type) const { return _ ## name [type]; }  \
  t& get ## name (AgentType type) { return _ ## name [type]; }  \
  void set ## name (AgentType type, const t& v) { _ ## name [type] = v; } \
  void set ## name (const t& v) { std::fill_n((_ ## name).begin(), (size_t)NAGENTS, v); }
#include "stat.def"

// --- Custom Accessors ---

  unsigned& getIntMtbFreq(size_t i) { return _intMtbFreq[i]; }
  unsigned getIntMtbFreq(size_t i) const { return _intMtbFreq[i]; }
  const std::vector<unsigned>& getIntMtbFreq() const { return _intMtbFreq; }
  void setIntMtbFreq(size_t i, const unsigned& v) { _intMtbFreq[i] = v; }

  Stat& getMacIntMtbStats(Mac::State i) { return _macIntMtbStats[i]; }
  const Stat& getMacIntMtbStats(Mac::State i) const { return _macIntMtbStats[i]; }

// --- Visitor Methods ---

  /*! \brief visit each stat
  * \tparam Visitor Visitor type that models the Visitor concept for each stat
  * \param v Visitor that visits each stat
  */
  template<typename Visitor>
  void visit(Visitor& v) {
    #define STAT(type, name, desc, reset)                v.visit(#name, _ ## name, desc);
    #define ARRAY_STAT(type, name, sz, desc, reset)      v.visit(#name, _ ## name, desc);
    #define AGENT_STAT(type, name, desc, reset)          ARRAY_STAT(type, name, NAGENTS, desc, reset)
    #define STATE_STAT(type, name, agent_t, desc, reset) ARRAY_STAT(type, name, agent_t::NSTATES, desc, reset)
    //#define AGENTSTATE_STAT(type, name, sz, desc, reset) v.visit(_ ## name); //???
    #include "stat.def"
  }

  /*! \brief visit each stat
  * \tparam Visitor Visitor type that models the Visitor concept for each stat
  * \param v Visitor that visits each stat
  */
  template<typename Visitor>
  void visit(Visitor& v) const {
    #define STAT(type, name, desc, reset)                v.visit(#name, _ ## name, desc);
    #define ARRAY_STAT(type, name, sz, desc, reset)      v.visit(#name, _ ## name, desc);
    #define AGENT_STAT(type, name, desc, reset)          ARRAY_STAT(type, name, NAGENTS, desc, reset)
    #define STATE_STAT(type, name, agent_t, desc, reset) ARRAY_STAT(type, name, agent_t::NSTATES, desc, reset)
    //#define AGENTSTATE_STAT(type, name, sz, desc, reset) v.visit(_ ## name); //???
    #include "stat.def"
  }

// --- Serialization ---

  /*! \brief serializes the stat object
  * \tparam Archive type that models boost::basearchive
  * \param ar Archive object to serialize to
  * \param file_version version of archive given
  */
  template<typename Archive>
  void serialize(Archive& ar, const unsigned int file_version = 0) {
    #define STAT(type, name, desc, reset)                ar & boost::serialization::make_nvp( #name, _ ## name);
    #define ARRAY_STAT(type, name, sz, desc, reset)      ar & boost::serialization::make_nvp( #name, _ ## name);
    #define AGENT_STAT(type, name, desc, reset)          ar & boost::serialization::make_nvp( #name, _ ## name);
    #define STATE_STAT(type, name, agent_t, desc, reset) ar & boost::serialization::make_nvp( #name, _ ## name);
    //#define AGENTSTATE_STAT(type, name, sz, desc, reset) ar & boost::serialization::make_nvp( #name, _ ## name); 
    #include "stat.def"
    ar & _intMtbFreq.size();
    for (size_t i = 0; i < _intMtbFreq.size(); i++)
      ar & _intMtbFreq[i];
  }

  void serialize(std::ostream& s, const unsigned int file_version = 0) const {
    #define STAT(type, name, desc, reset)           s << ( _ ## name) << std::endl;
    #define ARRAY_STAT(type, name, sz, desc, reset) s << ( _ ## name) << std::endl;
    #define AGENT_STAT(type, name, desc, reset)     s << ( _ ## name) << std::endl;
    #define STATE_STAT(type, name, agent_t, desc, reset) s << ( _ ## name) << std::endl;
    //#define AGENTSTATE_STAT(type, name, sz, desc, reset) s << boost::serialization::make_nvp( #name, _ ## name); 
    #include "stat.def"
    s << _intMtbFreq.size() << std::endl;
    for (size_t i = 0; i < _intMtbFreq.size(); i++)
      s << _intMtbFreq[i] << std::endl;
  }

  void clear() {
    #define STAT(type, name, desc, reset)           _##name = type();
    #define ARRAY_STAT(type, name, sz, desc, reset) std::fill_n(_##name.begin(), (size_t)sz, type());
    #define AGENT_STAT(type, name, desc, reset)     ARRAY_STAT(type, name, NAGENTS, desc, reset)
    #define STATE_STAT(type, name, agent_t, desc, reset) ARRAY_STAT(type, name, agent_t::NSTATES, desc, reset)
    //#define AGENTSTATE_STAT(type, name, sz, desc, reset) s << boost::serialization::make_nvp( #name, _ ## name); 
    #include "stat.def"
  }

  void reset() {
    #define STAT(type, name, desc, reset)           BOOST_PP_EXPR_IF(BOOST_PP_EQUAL(reset, 1), _##name = type());
    #define ARRAY_STAT(type, name, sz, desc, reset) BOOST_PP_EXPR_IF(BOOST_PP_EQUAL(reset, 1), std::fill_n(_##name.begin(), (size_t)sz, type()));
    #define AGENT_STAT(type, name, desc, reset)     ARRAY_STAT(type, name, NAGENTS, desc, reset)
    #define STATE_STAT(type, name, agent_t, desc, reset) ARRAY_STAT(type, name, agent_t::NSTATES, desc, reset)
    //#define AGENTSTATE_STAT(type, name, sz, desc, reset) s << boost::serialization::make_nvp( #name, _ ## name); 
    #include "stat.def"
  }
  void resetAgentStats() {
    #define STAT(type, name, desc, reset)           BOOST_PP_EXPR_IF(BOOST_PP_EQUAL(reset, 2), _##name = type());
    #define ARRAY_STAT(type, name, sz, desc, reset) BOOST_PP_EXPR_IF(BOOST_PP_EQUAL(reset, 2), std::fill_n(_##name.begin(), (size_t)sz, type()));
    #define AGENT_STAT(type, name, desc, reset)     ARRAY_STAT(type, name, NAGENTS, desc, reset)
    #define STATE_STAT(type, name, agent_t, desc, reset) ARRAY_STAT(type, name, agent_t::NSTATES, desc, reset)
    //#define AGENTSTATE_STAT(type, name, sz, desc, reset) s << boost::serialization::make_nvp( #name, _ ## name); 
    #include "stat.def"
    std::fill(_intMtbFreq.begin(), _intMtbFreq.end(), 0);
    std::fill_n(_macIntMtbStats.begin(), (size_t)Mac::NSTATES, Stat());
  }

  /*** Custom Functions ***/
  void updateAgentStatistics(const Agent* a) {
    a->updateStatistics(*this);
  }

};

template<>
inline void Stats::serialize(std::ostream& s, const unsigned int /* file_version */) {
  #define STAT(type, name, desc, reset)           s << ( _ ## name) << std::endl;
  #define ARRAY_STAT(type, name, sz, desc, reset) s << ( _ ## name) << std::endl;
  #define AGENT_STAT(type, name, desc, reset)     s << ( _ ## name) << std::endl;
  #define STATE_STAT(type, name, agent_t, desc, reset) s << ( _ ## name) << std::endl;
  //#define AGENTSTATE_STAT(type, name, sz, desc, reset) s << boost::serialization::make_nvp( #name, _ ## name); 
  #include "stat.def"
  s << _intMtbFreq.size() << std::endl;
  for (size_t i = 0; i < _intMtbFreq.size(); i++)
    s << _intMtbFreq[i] << std::endl;
}

template<>
inline void Stats::serialize(std::istream& s, const unsigned int /* file_version*/) {
  #define STAT(type, name, desc, reset)           s >> ( _ ## name);
  #define ARRAY_STAT(type, name, sz, desc, reset) s >> ( _ ## name);
  #define AGENT_STAT(type, name, desc, reset)     s >> ( _ ## name);
  #define STATE_STAT(type, name, agent_t, desc, reset) s >> ( _ ## name);
  //#define AGENTSTATE_STAT(type, name, sz, desc, reset) ar >> boost::serialization::make_nvp( #name, _ ## name); 
  #include "stat.def"
  unsigned tmp = 0;
  s >> tmp;
  assert(tmp == _intMtbFreq.size());
  for (size_t i = 0; i < _intMtbFreq.size(); i++)
    s >> _intMtbFreq[i];
}

#endif
