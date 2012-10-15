/*
 * tregulatory.h
 *
 *  Created on: 05-nov-2009
 *      Author: M. El-Kebir
 */

#ifndef TREGULATORY_H
#define TREGULATORY_H

#include "gr.h"
#include "tcell.h"

/**
* @brief
* @details Tregs look in their Moore neighborhood and deactivate all cells
* found with a probability (function of TNF and IL10 concentrations in the
* microcompartment).
* - Mr-: ALL OFF and deactivation time (no interaction with bug).
* - Mi-: ALL OFF and NO deactivation time (NF&kappa;B gives less TNF synthesis)
*   Bi still proliferates and Be still uptaken
* - Ma-: ALL OFF and deactivation time (no interaction with bug). Then the Ma
*   becomes Mr or Mi, with lifespan back to initial Mr.
* - Mci cannot be deactivated.
* - A deactivated cell can be re-activated again by resetting the clock of
*   deactivation time
*/
class Treg : public Tcell
{
public:
  enum State {TREG_ACTIVE, TREG_DEAD, NSTATES};
private:

  /*
   * !!! If the data members change then the serialize and deserialize functions need to be updated !!!
   */
  Treg::State _state;
  Treg::State _nextState;
  void handleResting(const int time, GrGrid& grid, Stats& stats);

protected:
  friend class boost::serialization::access;
  Treg();
public:
  Treg(int birthtime, int row, int col, Treg::State state);
  ~Treg();
  void move(GrGrid& grid);
  void secrete(GrGrid& grid, bool tnfrDynamics, bool nfkbDynamics, bool tnfDepletion, bool il10rDynamics, bool il10Depletion, double mdt);
  /**
  * @copydoc Agent::computeNextState
  * @details Tregs have only two explicit states: Active and Dead
  * \ref Treg::State.  Tregs can die via the following two mechanisms:
  * Apoptosis
  * :    TNF-dependent (through TNFR1 pathway, same as TNF apoptosis for macs)
  * Age
  * :    (3 days)
  */
  void computeNextState(const int time, GrGrid& grid, Stats& stats, bool tnfrDynamics, bool nfkbDynamics, bool il10Depletion, bool);
  void updateState();
  void updateStatistics(Stats& s) const;
  int getState() const;
  Treg::State getNextState() const;
  void kill();
  void deactivate(const int time, Stats& stats);
  bool isDead() const;
  bool isDeadNext();
  static bool isTreg(const Agent* pAgent);
  static bool isTreg(const Agent* pAgent, Treg::State state);
  void print() const;
  /*virtual*/
  Agent* clone() const
  {
    return new Treg(*this);
  }
  /**
  * @copydoc GrSimulation::serialize
  */
  template<class Archive>
  void serialize(Archive& ar, const unsigned int version);
  AgentType getAgentType() const;
};

inline AgentType Treg::getAgentType() const
{
  return TREG;
}

inline int Treg::getState() const
{
  return _state;
}

inline Treg::State Treg::getNextState() const
{
  return _nextState;
}

inline bool Treg::isTreg(const Agent* pAgent)
{
  return pAgent && pAgent->getAgentType() == TREG;
}

inline bool Treg::isTreg(const Agent* pAgent, Treg::State state)
{
  if (!isTreg(pAgent))
    {
      return false;
    }
  else
    {
      const Treg* pTreg = static_cast<const Treg*>(pAgent);
      return pTreg->getState() == state;
    }
}

inline bool Treg::isDead() const
{
  return _state == TREG_DEAD;
}

inline bool Treg::isDeadNext()
{
  return _nextState == TREG_DEAD;
}

std::ostream& operator<<(std::ostream& os, const Treg::State& s);

template<class Archive>
void Treg::serialize(Archive& ar, const unsigned int /*version*/) {
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Tcell);
  ar & boost::serialization::make_nvp("state", _state);
  ar & boost::serialization::make_nvp("nextState", _nextState);
}
#endif /* TREGULATORY_H */
