/*
 * tcytotoxic.h
 *
 *  Created on: 05-nov-2009
 *      Author: M. El-Kebir
 */

#ifndef TCYTOTOXIC_H
#define TCYTOTOXIC_H

#include "gr.h"
#include "tcell.h"

class Tcyt : public Tcell
{
public:
  enum State {TCYT_ACTIVE, TCYT_DOWN_REGULATED, TCYT_DEAD, NSTATES};
private:
  static const std::string _ClassName;

  /*
   * !!! If the data members change then the serialize and deserialize functions need to be updated !!!
   */
  Tcyt::State _state;
  Tcyt::State _nextState;
  int _deactivationTime;

  void handleActive(const int time, GrGrid& grid, Stats& stats);
  void handleDownRegulated(const int time, GrGrid& grid, Stats& stats);

protected:
  friend class boost::serialization::access;
  Tcyt();
public:
  Tcyt(int birthtime, int row, int col, Tcyt::State state);
  ~Tcyt();
  void move(GrGrid& grid);
  void secrete(GrGrid& grid, bool tnfrDynamics, bool nfkbDynamics, bool tnfDepletion, bool il10rDynamics, bool il10Depletion, double mdt);
  void computeNextState(const int time, GrGrid& grid, Stats& stats, bool tnfrDynamics, bool nfkbDynamics, bool il10rDynamics, bool);
  void updateState();
  void updateStatistics(Stats& s) const;
  int getState() const;
  Tcyt::State getNextState() const;
  void deactivate(const int time, Stats& stats);
  void kill();
  bool isDead() const;
  bool isDeadNext();
  int getDeactivationTime() const;
  static bool isTcyt(const Agent* pAgent);
  static bool isTcyt(const Agent* pAgent, Tcyt::State state);
  void print() const;
  /*virtual*/
  Agent* clone() const
  {
    return new Tcyt(*this);
  }
  /**
  * @copydoc GrSimulation::serialize
  */
  template<class Archive>
  void serialize(Archive& ar, const unsigned int version);
  AgentType getAgentType() const;
  /**
  * @copydoc Agent::visitProperties(Visitor& v)
  */
  template<typename Visitor>
  void visitProperties(Visitor& v);
  /**
  * @copydoc Agent::visitProperties(Visitor& v) const
  */
  template<typename Visitor>
  void visitProperties(Visitor& v) const;
};

inline AgentType Tcyt::getAgentType() const
{
  return TCYT;
}

inline int Tcyt::getDeactivationTime() const
{
  return _deactivationTime;
}

inline int Tcyt::getState() const
{
  return _state;
}

inline Tcyt::State Tcyt::getNextState() const
{
  return _nextState;
}

inline bool Tcyt::isTcyt(const Agent* pAgent)
{
  return pAgent && pAgent->getAgentType() == TCYT;
}

inline bool Tcyt::isTcyt(const Agent* pAgent, Tcyt::State state)
{
  if (!isTcyt(pAgent))
    {
      return false;
    }
  else
    {
      const Tcyt* pTcyt = static_cast<const Tcyt*>(pAgent);
      return pTcyt->getState() == state;
    }
}

inline bool Tcyt::isDead() const
{
  return _state == TCYT_DEAD;
}

inline bool Tcyt::isDeadNext()
{
  return _nextState == TCYT_DEAD;
}

template<typename Visitor>
inline void Tcyt::visitProperties(Visitor& v)
{
  v.visit("deactivationTime", _deactivationTime, "");
  Agent::visitProperties(v);
}
template<typename Visitor>
inline void Tcyt::visitProperties(Visitor& v) const
{
  v.visit("deactivationTime", _deactivationTime, "");
  Agent::visitProperties(v);
}

std::ostream& operator<<(std::ostream& os, const Tcyt::State& s);

template<class Archive>
void Tcyt::serialize(Archive& ar, const unsigned int /*version*/) {
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Tcell);
  ar & boost::serialization::make_nvp("state", _state);
  ar & boost::serialization::make_nvp("nextState", _nextState);
  ar & boost::serialization::make_nvp("deactivationTime", _deactivationTime);
}

#endif /* TCYTOTOXIC_H */
