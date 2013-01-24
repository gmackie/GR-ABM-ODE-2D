/*
 * tgamma.cpp
 *
 *  Created on: 03-nov-2009
 *      Author: M. El-Kebir
 */

#include "tgamma.h"
#include "macrophage.h"
#include "grgrid.h"
#include "stat.h"

using namespace std;

// Needed for deserializing the model state.
// Avoids the calls to the random number generator in the normal constructor, allowing the random number generator
// to remain in synch after deserialization.
Tgam::Tgam()
  : Tcell()
  , _state(TGAM_DEAD)
  , _nextState(TGAM_DEAD)
  , _deactivationTime(-1)
  , _transitionTime(-1)
  , _nAntigenStim(-1)
  , _nDownRegulated(-1)
  , _nICOS(-1)
{
}

Tgam::Tgam(int birthtime, int row, int col, Tgam::State state)

  : Tcell(birthtime, row, col, _PARAM(_kSynthTcell))
  , _state(state)
  , _nextState(state)
  , _deactivationTime(-1)
  , _transitionTime(-1)

  , _nAntigenStim(0)
  , _nDownRegulated(0)
  , _nICOS(0)
{
}

Tgam::~Tgam()
{
}

void Tgam::move(GrGrid& grid)
{
  Tcell::moveTcell(grid, true, true, true);
}

void Tgam::checkMacTgamInduced(GrGrid& grid)
{
    Pos coord;
    if (returnRandMacFromMoore(grid, coord))
    {
        Mac* pMac = dynamic_cast<Mac*>(grid.agent(coord, 0));
        if (!pMac) pMac = dynamic_cast<Mac*>(grid.agent(coord, 1));
        assert(pMac);

        // If the mac died on this time step ignore it.
        if (pMac->getNextState() == Mac::MAC_DEAD)
            return;

        // If mac does not die then check for induction conditions - Macs

        bool ProbAntigenPres = g_Rand.getReal() < _PARAM(Tcell_Tgam_probAntigenPresentation);

        if(pMac && pMac->getICOS())
            _nICOS ++;

        switch(pMac->getState())
         {
            case Mac::MAC_INFECTED:
            case Mac::MAC_CINFECTED:
            case Mac::MAC_ACTIVE:
                if (ProbAntigenPres)
                    _nAntigenStim ++;
                break;
            case Mac::MAC_RESTING:
            if ((pMac->getNFkB() || pMac->getStat1()) && ProbAntigenPres)
                    _nAntigenStim ++;
                break;
            throw std::runtime_error("Unknown Mac State in Tgam Induction Switch");
         }
    }
}

bool Tgam::checkTgamTransition(GrGrid& grid)
{
    checkMacTgamInduced(grid);

    double ProbSum = 0.0; // Initialize the probability sum for transition to TGAM_ACTIVE_DOUBLE

    if (_nAntigenStim > 1)  {
        double AgProb = (1.0/3.0) * (1 - exp((-_PARAM(Tcell_Tgam_rateAgDegree)*(_nAntigenStim - 1.0))));
        ProbSum += AgProb;
      }

    if (_nDownRegulated > 0)  {
        double TGFProb = (1.0/3.0) * (1 - exp((-_PARAM(Tcell_Tgam_rateTGFB)*(_nDownRegulated))));
        ProbSum += TGFProb;
      }

    if (_nICOS > 0)  {
        double ICOSProb = (1.0/3.0) * (1 - exp((-_PARAM(Tcell_Tgam_rateICOS)*(_nICOS))));
        ProbSum += ICOSProb;
      }

    assert(ProbSum >= 0 || ProbSum <= 1);
    return (g_Rand.getReal() < ProbSum);
}


void Tgam::secrete(GrGrid& grid, bool tnfrDynamics, bool, bool tnfDepletion, bool il10rDynamics, bool il10Depletion, double mdt)
{

    // This may not be necessary anymore since we can just use the TGAM_DOWN_REGULATED state in the switch
    if (_deactivationTime != -1)  {
      secTNF(grid, mdt, 0.0, _pos, tnfrDynamics, tnfDepletion, 0.0, 0.0, 0.0, _kmRNA, _kSynth);
      secIL10(grid, mdt, 0.0, _pos, il10rDynamics, il10Depletion, 0.0, 0.0, _kISynth);
      return;
    }

#if 0
  _kSynth = _PARAM(_kSynthTcell);
  _kmRNA = _PARAM(_kRNATcell);
//    calcIkmRNA(grid, _kmRNA, _kSynth, il10rDynamics);
#endif

  switch(_state)  {
  case Tgam::TGAM_ACTIVE_DOUBLE:
      secIL10(grid, mdt, 0.5, _pos, il10rDynamics, il10Depletion, _PARAM(_IkSynthTcell), (_PARAM(Tcell_Treg_dIL10)), _kISynth);
      secTNF(grid, mdt, 1.0, _pos, tnfrDynamics, tnfDepletion, _PARAM(_kSynthTcell), _PARAM(_kRNATcell), _PARAM(Tcell_Tgam_dTNF), _kmRNA, _kSynth);
      break;
  case Tgam::TGAM_INDUCED_REG:
      secIL10(grid, mdt, 1.0, _pos, il10rDynamics, il10Depletion, _PARAM(_IkSynthTcell), (_PARAM(Tcell_Treg_dIL10)), _kISynth);
      secTNF(grid, mdt, 1.0, _pos, tnfrDynamics, tnfDepletion, _PARAM(_kSynthTcell), _PARAM(_kRNATcell), _PARAM(Tcell_Tgam_dTNF), _kmRNA, _kSynth);
      break;
  case Tgam::TGAM_ACTIVE:
      secIL10(grid, mdt, 0.0, _pos, il10rDynamics, il10Depletion, _PARAM(_IkSynthTcell), _PARAM(Tcell_Tgam_dIL10), _kISynth);
      secTNF(grid, mdt, 1.0, _pos, tnfrDynamics, tnfDepletion, _PARAM(_kSynthTcell), _PARAM(_kRNATcell), _PARAM(Tcell_Tgam_dTNF), _kmRNA, _kSynth);
      break;
  case Tgam::TGAM_DOWN_REGULATED:
      secTNF(grid, mdt, 0.0, _pos, tnfrDynamics, tnfDepletion, 0.0, 0.0, 0.0, _kmRNA, _kSynth);
      secIL10(grid, mdt, 0.0, _pos, il10rDynamics, il10Depletion, 0.0, 0.0, _kISynth);
      break;
  case Tgam::TGAM_DEAD:
      break;
  default:
      throw std::runtime_error("Secretion Function - Invalid Tgam State");

  }
#if 0
  if (_state == TGAM_ACTIVE_DOUBLE)
    {
      _kISynth = 0.5 * _PARAM(_IkSynthTcell);

      if (!il10rDynamics && !il10Depletion)
        {
          grid.incil10(_pos, (0.5 * _PARAM(_dIL10_Treg) * mdt));
        }
    }

  else if (_state == TGAM_INDUCED_REG)
    {
      _kISynth = _PARAM(_IkSynthTcell);
      if (!il10rDynamics && !il10Depletion)
        {
          grid.incil10(_pos, (_PARAM(_dIL10_Treg) * mdt));
        }
    }

  else
    {
      _kISynth = 0;
    }

  if (!tnfrDynamics && !tnfDepletion)
    {
      double il10 = log(((grid.il10(_pos) * MW_IL10 * 1e6)/(NAV * VOL))); // converting il10 concentration to log(ng/mL) for use in dose dependence
      double tnfMOD = (1.0/(1.0 + exp((il10 + _PARAM(_LinkLogAlpha))/_PARAM(_LinkLogBeta)))); // calculate the fraction of inhibition
      grid.incTNF(_pos, (tnfMOD * _PARAM(_dTNF_Tgam) * mdt));
    }
#endif
}

void Tgam::computeNextState(const int time, GrGrid& grid, Stats& stats, bool tnfrDynamics, bool, bool, bool tgammatransition)
{
  // check if it is time to die
  if (timeToDie(time))  {
      _nextState = TGAM_DEAD;
    }

  // Always pass in false for nfkbDynamics for T cell apoptosis since they DO NOT have NFkB dynamics
  else if (TNFinducedApoptosis(grid, tnfrDynamics, false))  {
      ++stats.getTcellApoptosisTNF();
      _nextState = TGAM_DEAD;
      grid.incKillings(_pos);
    }
  else  {
      switch (_state)  {
        case TGAM_DEAD:
          // if dead, stay dead
          _nextState = TGAM_DEAD;
          break;
        case TGAM_ACTIVE:
          handleActive(time, grid, stats, tgammatransition);
          break;
        case TGAM_DOWN_REGULATED:
          handleDownRegulated(time, grid, stats);
          break;
        case TGAM_ACTIVE_DOUBLE:
          handleActiveDouble(time, grid, stats);
          break;
        case TGAM_INDUCED_REG:
          handleInducedReg(time, grid, stats);
          break;

        default:
          throw std::runtime_error("Unknown Tgam state");
        }
    }
}

void Tgam::handleActive(const int time, GrGrid& grid, Stats& stats, bool tgammatransition)
{
    if (g_Rand.getReal() < _PARAM(Tcell_Tgam_probApoptosisFasFasL))  {
        fasLigandKilling(grid, stats);
    }

    if(tgammatransition && checkTgamTransition(grid))  {
        _nextState = TGAM_ACTIVE_DOUBLE;
        _transitionTime = time;
     }
    else
        _nextState = TGAM_ACTIVE;
}

void Tgam::handleDownRegulated(const int time, GrGrid&, Stats&)
{
  if (time - _deactivationTime >= _PARAM(Tcell_Tgam_maxTimeReg))  {
      _nextState = TGAM_ACTIVE;
      _deactivationTime = -1;
    }
  else  {
      _nextState = TGAM_DOWN_REGULATED;
    }
}

void Tgam::handleActiveDouble(const int time, GrGrid& grid, Stats& stats)
{
    // Carries out same action as TGAM_ACTIVE but secretes half rate of IL10
    // Distinct state so it is easier to track

    if (time - _transitionTime >= _PARAM(Tcell_Tgam_maxTimeDouble))  {
        _nextState = TGAM_INDUCED_REG;
      }
    else  {
        _nextState = TGAM_ACTIVE_DOUBLE;
      }

    if (g_Rand.getReal() < _PARAM(Tcell_Tgam_probApoptosisFasFasL))  {
        fasLigandKilling(grid, stats);
    }
}

void Tgam::handleInducedReg(const int, GrGrid&, Stats&)
{
  _nextState = TGAM_INDUCED_REG;
}

void Tgam::fasLigandKilling(GrGrid& grid, Stats& stats)
{
    Pos coord;
    if (returnRandMacFromMoore(grid, coord))  {
        Mac* pMac = dynamic_cast<Mac*>(grid.agent(coord, 0));
        if (!pMac) pMac = dynamic_cast<Mac*>(grid.agent(coord, 1));
        assert(pMac);

        // If the mac died on this time step ignore it.
        if (pMac->getNextState() == Mac::MAC_DEAD)  {
            return;
          }

        // Fas/FasL induced apoptosis with probability
        if (pMac->getState() == Mac::MAC_INFECTED || pMac->getState() == Mac::MAC_CINFECTED)  {
            ++stats.getApoptosisFasFasL();
            pMac->apoptosis(grid);
            pMac->kill();
            grid.incKillings(coord);
          }
    }
}

void Tgam::deactivate(const int time, Stats& stats)
{

  switch (_state)  {
    case Tgam::TGAM_INDUCED_REG:
      break;
    case Tgam::TGAM_DOWN_REGULATED:
      break;
    case Tgam::TGAM_ACTIVE:
    case Tgam::TGAM_ACTIVE_DOUBLE:
      ++stats.getTgamDeactivation(_state);
      _nDownRegulated ++;
      _nextState = _state = TGAM_DOWN_REGULATED;
      _deactivationTime = time;
      break;
    default:
      ;
    }
}

void Tgam::updateState()
{
  _state = _nextState;
}

void Tgam::kill()
{
  _nextState = _state = TGAM_DEAD;
}

void Tgam::print() const
{
  std::cout << "Tgam - " << _birthTime << " - ";

  switch (_state)
    {
    case TGAM_DEAD:
      std::cout << "dead" << std::endl;
      break;
    case TGAM_DOWN_REGULATED:
      std::cout << "down-regulated" << std::endl;
      break;
    case TGAM_ACTIVE:
      std::cout << "active" << std::endl;
      break;
    default:
      throw std::runtime_error("Unknown Tgam state");
    }
}

#if 0

void Tgam::serialize(std::ostream& out) const
{
  assert(out.good());

  Serialization::writeHeader(out, Tgam::_ClassName);

  Tcell::serialize(out);

  int intVal = (int) _state;
  out << intVal << std::endl;

  intVal = (int) _nextState;
  out << intVal << std::endl;

  out << _deactivationTime << std::endl;
  out << _transitionTime << std::endl;

  out << _nAntigenStim << std::endl;
  out << _nDownRegulated << std::endl;
  out << _nICOS << std::endl;

  Serialization::writeFooter(out, Tgam::_ClassName);
}

void Tgam::deserialize(std::istream& in)
{
  assert(in.good());

  if (!Serialization::readHeader(in, Tgam::_ClassName))
    {
      exit(1);
    }

  int intVal;

  Tcell::deserialize(in);

  in >> intVal;
  _state = (Tgam::State) intVal;

  in >> intVal;
  _nextState = (Tgam::State) intVal;

  in >> _deactivationTime;
  in >> _transitionTime;

  in >> _nAntigenStim;
  in >> _nDownRegulated;
  in >> _nICOS;

  if (!Serialization::readFooter(in, Tgam::_ClassName))
    {
      exit(1);
    }
}

#endif

void Tgam::updateStatistics(Stats& s) const
{
  ++s.getNrOfAgents(TGAM);
  ++s.getNrOfTgams((Tgam::State)getState());
}

const char* _tgamStrings[] = { "Active", "Down Regulated", "Active Double", "Induced Regulated", "Dead"};

std::ostream& operator<<(std::ostream& os, const Tgam::State& s)
{
  return os << _tgamStrings[s];
}
