/*
 * tregulatory.cpp
 *
 *  Created on: 05-nov-2009
 *      Author: M. El-Kebir
 */

#include "tregulatory.h"
#include "grgrid.h"
#include "stat.h"

using namespace std;

// Needed for deserializing the model state.
// Avoids the calls to the random number generator in the normal constructor, allowing the random number generator
// to remain in synch after deserialization.
Treg::Treg()
  : Tcell()
  , _state(TREG_DEAD)
  , _nextState(TREG_DEAD)
{
}

Treg::Treg(int birthtime, int row, int col, Treg::State state)

  : Tcell(birthtime, row, col, 0.0 /* kSynth */)
  , _state(state)
  , _nextState(state)

{
}

Treg::~Treg()
{
}

void Treg::move(GrGrid& grid)
{
  Tcell::moveTcell(grid, false, true, false);
}

void Treg::secrete(GrGrid& grid, bool, bool, bool, bool il10rDynamics, bool il10Depletion, double mdt)
{

  secIL10(grid, mdt, 1.0, _pos, il10rDynamics, il10Depletion, _PARAM(_IkSynthTcell), (_PARAM(Tcell_Treg_dIL10)), _kISynth);

 #if 0
    _kISynth = _PARAM(_IkSynthTcell);

  if (!il10rDynamics && !il10Depletion)
    {
      grid.incil10(_pos, (_PARAM(_dIL10_Treg) * mdt));
    }
#endif

}

void Treg::deactivate(const int, Stats&)
{
}

void Treg::computeNextState(const int time, GrGrid& grid, Stats& stats, bool tnfrDynamics, bool, bool, bool)
{
  // check if it is time to die
  if (timeToDie(time))  {
      _nextState = TREG_DEAD;
    }

  // Always pass in false for nfkbDynamics for T cell apoptosis since they DO NOT have NFkB dynamics
  else if (TNFinducedApoptosis(grid, tnfrDynamics, false))  {
      ++stats.getTcellApoptosisTNF();
      _nextState = TREG_DEAD;
      grid.incKillings(_pos);
    }

  else  {
      switch (_state)
        {
        case TREG_DEAD:
          // if dead, stay dead
          _nextState = TREG_DEAD;
          break;
        case TREG_ACTIVE:
          handleResting(time, grid, stats);
          break;
        default:
          throw std::runtime_error("Unknown Tcyt state");
        }
    }
}

// Don't deactivate an agent that is already dead, either now or for the next time step -
// died of old age, apoptosis, etc.
// Otherwise its state and next state are set to deactivated, which essentially resurrects it.
// This is bad in general. It can also cause a problem if the resurrected agent is in a caseated
// compartment and the agent is serialized to a saved simulation state. If the saved state is
// loaded then the agent will be deserialized and added to a cell list, but it will fail to be
// added to its grid compartment because an agent can't be added to a caseated compartment.
// This will leave the simulation in an inconsistent state.

void Treg::handleResting(const int time, GrGrid& grid, Stats& stats)
{
#if 0
  Scalar IL10toTNFWeight = (_PARAM(_IkOn) * (_surfIL10R + _surfBoundIL10R)) / (_PARAM(_kOn1) * (_surfBoundTNFR1 + _surfTNFR1)); // Comparing Kd allows us to scale the bound receptors such that the numbers do not favor TNFs higher affinity
  Scalar numberFractionTNF;
  if (_surfBoundTNFR1 == 0.0 && _surfBoundIL10R == 0.0)
    numberFractionTNF = 0.0;
  else
    numberFractionTNF = (IL10toTNFWeight * _surfBoundTNFR1)/((IL10toTNFWeight * _surfBoundTNFR1) + _surfBoundIL10R); // Calculate the scaled number fraction of TNF (Like mol/weight fraction)

//    std::cout << "Pos: " << _pos << "  Xtnf: " << numberFractionTNF << std::endl;
#endif

  for (int i = -1; i <= 1; i++)
    {
      for (int j = -1; j <= 1; j++)
        {
          Pos p(grid.mod_row(_pos.x+i), grid.mod_col(_pos.y+j));

          for(unsigned k=0; k<GrGrid::MAX_AGENTS_PER_CELL; k++)
            {
              Agent* pAgent = grid.agent(p, k);
              if(pAgent && !pAgent->isDead() && !pAgent->isDeadNext())
                {
                  Scalar coinFlip = g_Rand.getReal();
                  if (calcDeactivationProbability(coinFlip, _surfIL10R, _surfBoundIL10R, _surfTNFR1, _surfBoundTNFR1))  {
                      grid.agent(p, k)->deactivate(time, stats);
                  }
#if 0
                  Scalar scaledProb = _PARAM(Tcell_Treg_probTregDeactivate) * ((numberFractionTNF * _PARAM(Tcell_Treg_deactivateSlope)) + _PARAM(Tcell_Treg_deactivateIntercept));
//                      std::cout << "Pos :" << _pos << "  IL10: " << grid.il10(_pos) << "  TNF: " << grid.TNF(_pos) << "  Scaled Prob: " << scaledProb << std::endl;
                  if (coinFlip  <= scaledProb)
                    {
                      grid.agent(p, k)->deactivate(time, stats);
                    }
#endif
                }
            }
        }
    }
}

bool Treg::calcDeactivationProbability(const Scalar rnum, const Scalar il10r, const Scalar boundil10r, const Scalar tnfr1, const Scalar boundtnfr1)
{

    // This function roughly simulates the deactivation mechanisms that are
    // not occuring through IL-10. We use boundTNFR1 and boundIL10R to scale the
    // probability of deactivation such that if there is a large amount of IL-10 present
    // the probability is lowered, while if a large amount of TNF is present the probability is raised.
    // This allows for 'compensation' of other deactivation mechanisms when IL-10 is not present.

    Scalar numberFractionTNF;

    // Calculate the weight of IL10 to TNF based on total receptor number
    // and the relative strength of binding. This prevents one receptor system
    // from dominating the calculation if it has more receptors or binds faster
    Scalar IL10toTNFWeight = (_PARAM(_IkOn) * (il10r + boundil10r)) / (_PARAM(_kOn1) * (tnfr1 + boundtnfr1));

    if (boundtnfr1 == 0.0 && boundil10r == 0.0)
      numberFractionTNF = 0.0;
    else
      numberFractionTNF = (IL10toTNFWeight * boundtnfr1)/((IL10toTNFWeight * boundtnfr1) + boundil10r); // Calculate the scaled number fraction of boundTNFR1 (Like mol/weight fraction)

    Scalar scaledProb = _PARAM(Tcell_Treg_probTregDeactivate) * ((numberFractionTNF * _PARAM(Tcell_Treg_deactivateSlope)) + _PARAM(Tcell_Treg_deactivateIntercept));

    return (rnum <= scaledProb);
}

void Treg::updateState()
{
  _state = _nextState;
}

void Treg::kill()
{
  _nextState = _state = TREG_DEAD;
}

void Treg::print() const
{
  std::cout << "Treg - " << _birthTime << " - ";

  switch (_state)
    {
    case TREG_DEAD:
      std::cout << "dead" << std::endl;
      break;
    case TREG_ACTIVE:
      std::cout << "active" << std::endl;
      break;
    default:
      throw std::runtime_error("Unknown Tcyt state");
    }
}

#if 0
void Treg::serialize(std::ostream& out) const
{
  assert(out.good());

  Serialization::writeHeader(out, Treg::_ClassName);

  Tcell::serialize(out);

  int intVal = (int) _state;
  out << intVal << std::endl;

  intVal = (int) _nextState;
  out << intVal << std::endl;

  Serialization::writeFooter(out, Treg::_ClassName);
}

void Treg::deserialize(std::istream& in)
{
  assert(in.good());

  if (!Serialization::readHeader(in, Treg::_ClassName))
    {
      exit(1);
    }

  int intVal;

  Tcell::deserialize(in);

  in >> intVal;
  _state = (Treg::State) intVal;

  in >> intVal;
  _nextState = (Treg::State) intVal;

  if (!Serialization::readFooter(in, Treg::_ClassName))
    {
      exit(1);
    }
}
#endif

void Treg::updateStatistics(Stats& s) const
{
  ++s.getNrOfAgents(TREG);
  ++s.getNrOfTregs((Treg::State)getState());
}

const char* _tregStrings[] = { "Active", "Dead" };

std::ostream& operator<<(std::ostream& os, const Treg::State& s)
{
  return os << _tregStrings[s];
}
