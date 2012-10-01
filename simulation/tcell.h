/*
 * tcell.h
 *
 *  Created on: 13-nov-2009
 *      Author: M. El-Kebir
 */

#ifndef TCELL_H
#define TCELL_H

#include "agent.h"

class Tcell : public Agent
{
private:

  static int _tcellodeSize;
  /*
   * !!! If the data members change then the serialize and deserialize functions need to be updated !!!
   */
public:
  Tcell();
  Tcell(int birthTime, int row, int col, Scalar kSynth);

  virtual ~Tcell();
  void moveTcell(GrGrid& grid, bool ccl2, bool ccl5, bool cxcl9);
  static void setTcellOdeSize(int odesize);
  static bool isTcell(const Agent* pAgent);
  template<class Archive>
  void serialize(Archive& ar, const unsigned int version);
  virtual void solveDegradation (GrGrid& grid, double dt, bool tnfrDynamics, bool il10rDynamics);

};

inline bool Tcell::isTcell(const Agent* pAgent)
{
  return pAgent && pAgent->getAgentType() != MAC;
}

inline void Tcell::setTcellOdeSize(int odesize)
{
  _tcellodeSize = odesize;
}

template<class Archive>
void Tcell::serialize(Archive& ar, const unsigned int /*version*/) {
  ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Agent);
  ar & BOOST_SERIALIZATION_NVP(_tcellodeSize);
}

BOOST_SERIALIZATION_ASSUME_ABSTRACT(Tcell);
#endif /* TCELL_H */
