/*
 * tregulatory.cpp
 *
 *  Created on: 05-nov-2009
 *      Author: M. El-Kebir
 */

#include "tregulatory.h"
#include "grgrid.h"
#include "stat.h"
#include "serialization.h"

using namespace std;

const string Treg::_ClassName = "Treg";


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

  _kISynth = _PARAM(PARAM_GR_I_K_SYNTH_TCELL);

  if (!il10rDynamics && !il10Depletion)
    {
      grid.incil10(_pos, (_PARAM(PARAM_TREG_SEC_RATE_IL10) * mdt));
    }

}

void Treg::deactivate(const int, Stats&)
{
}

void Treg::computeNextState(const int time, GrGrid& grid, Stats& stats, bool tnfrDynamics, bool, bool, bool)
{
  // check if it is time to die
  if (timeToDie(time))
    {
      _nextState = TREG_DEAD;
    }

  // Always pass in false for nfkbDynamics for T cell apoptosis since they DO NOT have NFkB dynamics
  else if (TNFinducedApoptosis(grid, tnfrDynamics, false))
    {
      ++stats.getTcellApoptosisTNF();
      _nextState = TREG_DEAD;
      grid.incKillings(_pos);
    }

//	else if (tnfrDynamics && intCompareGT(_intBoundTNFR1, _PARAM(PARAM_GR_THRESHOLD_APOPTOSIS_TNF_MOLECULAR)) &&
//			 intCompareGT(1 - pow(2.7183, -_PARAM(PARAM_GR_K_APOPTOSIS_MOLECULAR) * (_intBoundTNFR1 - _PARAM(PARAM_GR_THRESHOLD_APOPTOSIS_TNF_MOLECULAR))), g_Rand.getReal()))
//	{
//		// TNF induced apoptosis
//		++stats.getTcellApoptosisTNF();
//		_nextState = TREG_DEAD;
//        grid.incKillings(_pos);
//	}
//	else if (!tnfrDynamics && tnfBoundFraction > _PARAM(PARAM_GR_THRESHOLD_APOPTOSIS_TNF) &&
//			 g_Rand.getReal() < 1 - pow(2.7183, -_PARAM(PARAM_GR_K_APOPTOSIS) * (tnfBoundFraction - _PARAM(PARAM_GR_THRESHOLD_APOPTOSIS_TNF))))
//	{
//		// TNF induced apoptosis
//		++stats.getTcellApoptosisTNF();
//		_nextState = TREG_DEAD;
//        grid.incKillings(_pos);
//	}

  else
    {
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
//    Scalar currentTNF = grid.TNF(_pos); // Get TNF concentration for scaling Treg deactivation
//    Scalar currentIL10 = grid.il10(_pos); // Get IL10 concentration for scaling Treg deactivation

  Scalar IL10toTNFWeight = (_PARAM(PARAM_GR_I_K_ON) * (_surfIL10R + _surfBoundIL10R)) / (_PARAM(PARAM_GR_K_ON1) * (_surfBoundTNFR1 + _surfTNFR1)); // Comparing Kd allows us to scale the bound receptors such that the numbers do not favor TNFs higher affinity
  Scalar numberFractionTNF;
  if (_surfBoundTNFR1 == 0.0 && _surfBoundIL10R == 0.0)
    numberFractionTNF = 0.0;
  else
    numberFractionTNF = (IL10toTNFWeight * _surfBoundTNFR1)/((IL10toTNFWeight * _surfBoundTNFR1) + _surfBoundIL10R); // Calculate the scaled number fraction of TNF (Like mol/weight fraction)

//    std::cout << "Pos: " << _pos << "  Xtnf: " << numberFractionTNF << std::endl;

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

                  // Scale Treg deactivation by TNF and IL10 concentration (only monitor molecules we have) as a proxy for
                  // feedback to TGF-b cell-cell contact mechanism
                  // TNF is a MM type eqn
                  // IL10 is an exponential decay
                  // These are inverse equations of each other since TNF would upregulate the 'Monitor Molecule' while IL10 would
                  // downregulate the 'Monitor Molecule'
//                      Scalar scaledProbTNF = ((_PARAM(PARAM_TREG_PROB_DEACTIVATE) * currentTNF)/(currentTNF + _PARAM(PARAM_TREG_DEACTIVATE_HALF_SAT_TNF)));
//                      Scalar scaledProbIL10 = _PARAM(PARAM_TREG_PROB_DEACTIVATE) * pow(2.7183, -_PARAM(PARAM_TREG_DEACTIVATE_HALF_SAT_IL10) * currentIL10);

//                      std::cout << "Scaled Prob: " << scaledProbTNF + scaledProbIL10 << std::endl;

//                      Scalar scaledProb = std::min((scaledProbTNF + scaledProbIL10), _PARAM(PARAM_TREG_PROB_DEACTIVATE));

                  Scalar scaledProb = _PARAM(PARAM_TREG_PROB_DEACTIVATE) * ((numberFractionTNF * _PARAM(PARAM_TREG_DEACTIVATE_SLOPE)) + _PARAM(PARAM_TREG_DEACTIVATE_INTERCEPT));
//                      std::cout << "Pos :" << _pos << "  IL10: " << grid.il10(_pos) << "  TNF: " << grid.TNF(_pos) << "  Scaled Prob: " << scaledProb << std::endl;
                  if (coinFlip  <= scaledProb)
                    {
                      grid.agent(p, k)->deactivate(time, stats);
                    }
                }
            }
        }
    }
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
