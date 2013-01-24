/*
 * tcell.cpp
 *
 *  Created on: 13-nov-2009
 *      Author: M. El-Kebir
 */


#include "tcell.h"
#include "stat.h"
#include "grgrid.h"

using namespace std;

// Needed for deserializing the model state.
// Avoids the calls to the random number generator in the normal constructor, allowing the random number generator
// to remain in synch after deserialization.
Tcell::Tcell()
  : Agent()
{
}

Tcell::Tcell(int birthtime, int row, int col, Scalar kSynth)
  : Agent(birthtime, birthtime + _PARAM(Tcell_maxAge), row, col
          //TNFR Components
          , (Scalar) _PARAM(_meanTNFR1Tcell)
          , (Scalar) _PARAM(_stdTNFR1Tcell)
          , (Scalar) _PARAM(_meanTNFR2Tcell)
          , (Scalar) _PARAM(_stdTNFR2Tcell)
          , kSynth
          , (Scalar) _PARAM(_kTaceTcell)

          // IL10 components
          , (Scalar) _PARAM(_meanIL10RTcell)
          , (Scalar) _PARAM(_stdIL10RTcell)
         )
{
}

Tcell::~Tcell()
{
}

void Tcell::solveDegradation(GrGrid& grid, double dt, bool tnfrDynamics, bool il10rDynamics)
{
  //DBG
//	cout << "Debug: Running T cell degradation"
//		<<  " _PARAM(_meanTNFR1Tcell): " << _PARAM(_meanTNFR1Tcell)
//		<<  " _PARAM(_meanIL10RTcell): " << _PARAM(_meanIL10RTcell)
//		<< std::endl;
  // DBG

  Agent::solveDegradation(grid, dt, tnfrDynamics, il10rDynamics, _PARAM(_meanTNFR1Tcell), _PARAM(_meanIL10RTcell));
}

void Tcell::moveTcell(GrGrid& grid, bool ccl2, bool ccl5, bool cxcl9)
{
  Pos pos  = Agent::moveAgent(grid, ccl2, ccl5, cxcl9, false, _PARAM(Tcell_movementBonusFactor));

  // Check whether newCell is not caseated and contains empty slots
  size_t n = grid.getNumberOfAgents(pos);
  if (n != 2 && !grid.isCaseated(pos))
    {
      // Move with p = 1, if newCell is empty
      // Move with p = _PROB_MOVE_TCELL_TO_MAC, if newCell contains a macrophage
      // Move with p = _PROB_MOVE_TCELL_TO_TCELL, if newCell contains a T cell
      if ((n == 0) ||
          (grid.hasAgentType(MAC, pos) && g_Rand.getReal() < _PARAM(Tcell_probMoveToMac)) ||
          (grid.hasTcell(pos) && g_Rand.getReal() < _PARAM(Tcell_probMoveToTcell)))
        {
          assert_res(grid.removeAgent(this));
          grid.addAgent(this, pos);
          _pos = pos;
        }
    }
}

bool Tcell::returnRandMacFromMoore(GrGrid& grid, Pos& vectorPos)
{
    vector<int> PossibleOrdinal;

    for (int k=0; k<9; k++)  {
        Pos p(this->compartmentOrdinalToCoordinates(k, grid.getRange()));
        if(grid.hasAgentType(MAC, p))  {
            PossibleOrdinal.push_back(k);
          }
      }

    // If there are no Macs then do not kill anything
    if ((int) PossibleOrdinal.size() == 0)  {
        return false;
      }

    int PossibleOrds = PossibleOrdinal.size();
    int RandNum = g_Rand.getInt(PossibleOrds, 0);

    vectorPos = this->compartmentOrdinalToCoordinates(PossibleOrdinal[RandNum], grid.getRange());

    return true;
}
