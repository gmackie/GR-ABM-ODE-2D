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
                GR_DISSEMINATION, GR_DISSEMINATION_INCONSISTENT, GR_UNKNOWN, GR_NONE
             } GrStatus;
namespace ba = boost::accumulators;


inline std::ostream& operator<<(std::ostream& s, const GrStatus& a)
{
  s<<(int)a;
  return s;
}
inline std::istream& operator>>(std::istream& s, GrStatus& a)
{
  int x;
  s>>x;
  a = (GrStatus)x;
  return s;
}

// --- iostream members for boost::array<> ---
template<typename T, size_t N>
inline std::ostream& operator<<(std::ostream& s, const boost::array<T, N>& a)
{
  s<<N<<' ';
  for(typename boost::array<T, N>::const_iterator it = a.begin(); it != a.end(); it++)
    s<<*it<<' ';
  return s;
}

template<typename T, size_t N>
inline std::istream& operator>>(std::istream& s, boost::array<T, N>& a)
{
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
  typedef ba::accumulator_set<Scalar, ba::stats< ba::features< ba::tag::variance, ba::tag::sum, ba::tag::min, ba::tag::max, ba::tag::median > > > Stat;
  struct group_type {}; //For identifying groups
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
  Stat _macGrowthRateStat;
  std::vector<unsigned> _intMtbFreq;
  std::vector<unsigned> _growthRateFreq;

  // Total drug concentration on the grid is kept for a number of time steps
  // specified as a parameter file parameter, so a Boost array cannot be used
  // because that requires specifying the size at compile time.
  std::vector<Scalar> _drugConcentrationINH;

public:

// --- Constructor ---

  Stats() :
      _intMtbFreq(int(_PARAM(Mac_nrIntMtbBurstCInf)) + 1, 0)
    , _growthRateFreq(_PARAM(_growthRateSamples), 0)
    , _drugConcentrationINH(_PARAM(_drugConcentrationStatInterval), 0.0)
  {
    clear();
  }

  /// @name Accessors
  /// @{
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
  /// @}

// --- Custom Accessors ---

  Stat& getMacIntMtbStats(Mac::State i)
  {
    return _macIntMtbStats[i];
  }
  const Stat& getMacIntMtbStats(Mac::State i) const
  {
    return _macIntMtbStats[i];
  }

  Stat& getMacGrowthRateStat() {
    return _macGrowthRateStat;
  }
  const Stat& getMacGrowthRateStat() const {
    return _macGrowthRateStat;
  }

  unsigned& getIntMtbFreq(size_t i)
  {
    return _intMtbFreq[i];
  }
  unsigned getIntMtbFreq(size_t i) const
  {
    return _intMtbFreq[i];
  }
  const std::vector<unsigned>& getIntMtbFreq() const
  {
    return _intMtbFreq;
  }
  void setIntMtbFreq(size_t i, const unsigned& v)
  {
    _intMtbFreq[i] = v;
  }

  unsigned& getGrowthRateFreq(size_t i)
  {
    return _growthRateFreq[i];
  }
  unsigned getGrowthRateFreq(size_t i) const
  {
    return _growthRateFreq[i];
  }
  const std::vector<unsigned>& getGrowthRateFreq() const
  {
    return _growthRateFreq;
  }
  void setGrowthRateFreq(size_t i, const unsigned& v)
  {
    _growthRateFreq[i] = v;
  }

  Scalar& getDrugConcentrationINH(size_t i)
  {
    return _drugConcentrationINH[i];
  }

  Scalar getDrugConcentrationINH(size_t i) const
  {
    return _drugConcentrationINH[i];
  }

  const std::vector<Scalar>& getDrugConcentrationINH() const
  {
    return _drugConcentrationINH;
  }

  void setDrugConcentrationINH(size_t i, const Scalar v)
  {
    _drugConcentrationINH[i] = v;
  }



// --- Visitor Methods ---

  /*! \brief visit each stat
  * \tparam Visitor Visitor type that models the Visitor concept for each stat
  * \param v Visitor that visits each stat
  */
  template<typename Visitor>
  void visit(Visitor& v)
  {
      group_type g;
#define GROUP_STATS(name, desc) v.visit(name, g, desc);
#define STAT(type, name, desc, reset)                v.visit(#name, _ ## name, desc);
#define ARRAY_STAT(type, name, sz, desc, reset)      v.visit(#name, _ ## name, desc, sz);
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
  void visit(Visitor& v) const
  {
      group_type g;
#define GROUP_STATS(name, desc) v.visit(name, g, desc);
#define STAT(type, name, desc, reset)                v.visit(#name, _ ## name, desc);
#define ARRAY_STAT(type, name, sz, desc, reset)      v.visit(#name, _ ## name, desc, sz);
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
  void serialize(Archive& ar, const unsigned int file_version=0)
  {
    (void)file_version;
#define STAT(type, name, desc, reset)                ar & boost::serialization::make_nvp( #name, _ ## name);
#define ARRAY_STAT(type, name, sz, desc, reset)      ar & boost::serialization::make_nvp( #name, _ ## name);
#define AGENT_STAT(type, name, desc, reset)          ar & boost::serialization::make_nvp( #name, _ ## name);
#define STATE_STAT(type, name, agent_t, desc, reset) ar & boost::serialization::make_nvp( #name, _ ## name);
    //#define AGENTSTATE_STAT(type, name, sz, desc, reset) ar & boost::serialization::make_nvp( #name, _ ## name);
#include "stat.def"
    ar & boost::serialization::make_nvp("intMtbFreq", _intMtbFreq);
  }

  void serialize(std::ostream& s) const
  {
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

  /**
  * @brief Clears all stats with the value of the default constructor
  */
  void clear()
  {
#define STAT(type, name, desc, reset)           _##name = type();
#define ARRAY_STAT(type, name, sz, desc, reset) std::fill_n(_##name.begin(), (size_t)sz, type());
#define AGENT_STAT(type, name, desc, reset)     ARRAY_STAT(type, name, NAGENTS, desc, reset)
#define STATE_STAT(type, name, agent_t, desc, reset) ARRAY_STAT(type, name, agent_t::NSTATES, desc, reset)
    //#define AGENTSTATE_STAT(type, name, sz, desc, reset) s << boost::serialization::make_nvp( #name, _ ## name);
#include "stat.def"
  }

  /**
  * @brief Resets all resetable stats (reset=1) to the value of the default constructor
  */
  void reset()
  {
#define STAT(type, name, desc, reset)           BOOST_PP_EXPR_IF(BOOST_PP_EQUAL(reset, 1), _##name = type());
#define ARRAY_STAT(type, name, sz, desc, reset) BOOST_PP_EXPR_IF(BOOST_PP_EQUAL(reset, 1), std::fill_n(_##name.begin(), (size_t)sz, type()));
#define AGENT_STAT(type, name, desc, reset)     ARRAY_STAT(type, name, NAGENTS, desc, reset)
#define STATE_STAT(type, name, agent_t, desc, reset) ARRAY_STAT(type, name, agent_t::NSTATES, desc, reset)
    //#define AGENTSTATE_STAT(type, name, sz, desc, reset) s << boost::serialization::make_nvp( #name, _ ## name);
#include "stat.def"
  }
  /**
  * @brief Resets all resetable stats (reset=2) to the value of the default constructor
  */
  void resetAgentStats()
  {
#define STAT(type, name, desc, reset)           BOOST_PP_EXPR_IF(BOOST_PP_EQUAL(reset, 2), _##name = type());
#define ARRAY_STAT(type, name, sz, desc, reset) BOOST_PP_EXPR_IF(BOOST_PP_EQUAL(reset, 2), std::fill_n(_##name.begin(), (size_t)sz, type()));
#define AGENT_STAT(type, name, desc, reset)     ARRAY_STAT(type, name, NAGENTS, desc, reset)
#define STATE_STAT(type, name, agent_t, desc, reset) ARRAY_STAT(type, name, agent_t::NSTATES, desc, reset)
    //#define AGENTSTATE_STAT(type, name, sz, desc, reset) s << boost::serialization::make_nvp( #name, _ ## name);
#include "stat.def"
    std::fill_n(_macIntMtbStats.begin(), (size_t)Mac::NSTATES, Stat());
    _macGrowthRateStat = Stat();
    std::fill(_intMtbFreq.begin(), _intMtbFreq.end(), 0);
    std::fill(_growthRateFreq.begin(), _growthRateFreq.end(), 0);
  }

  /*** Custom Functions ***/
  /**
  * @brief Double dispatch method for updating agent statistics
  *
  * @param a
  */
  void updateAgentStatistics(const Agent* a)
  {
    a->updateStatistics(*this);
  }

};

template<>
inline void Stats::serialize(std::ostream& s, const unsigned int /* file_version */)
{
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
inline void Stats::serialize(std::istream& s, const unsigned int /* file_version*/)
{
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
