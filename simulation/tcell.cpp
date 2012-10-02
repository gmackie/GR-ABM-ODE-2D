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
  : Agent(birthtime, birthtime + _PARAM(PARAM_TCELL_AGE), row, col
          //TNFR Components
          , (Scalar) _PARAM(PARAM_GR_MEAN_TNFR1_TCELL)
          , (Scalar) _PARAM(PARAM_GR_STD_TNFR1_TCELL)
          , (Scalar) _PARAM(PARAM_GR_MEAN_TNFR2_TCELL)
          , (Scalar) _PARAM(PARAM_GR_STD_TNFR2_TCELL)
          , kSynth
          , (Scalar) _PARAM(PARAM_GR_K_TACE_TCELL)

          // IL10 components
          , (Scalar) _PARAM(PARAM_GR_I_IL10R_TCELL)
          , (Scalar) _PARAM(PARAM_GR_STD_IL10R_TCELL)
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
//		<<  " _PARAM(PARAM_GR_MEAN_TNFR1_TCELL): " << _PARAM(PARAM_GR_MEAN_TNFR1_TCELL)
//		<<  " _PARAM(PARAM_GR_I_IL10R_TCELL): " << _PARAM(PARAM_GR_I_IL10R_TCELL)
//		<< std::endl;
  // DBG

  Agent::solveDegradation(grid, dt, tnfrDynamics, il10rDynamics, _PARAM(PARAM_GR_MEAN_TNFR1_TCELL), _PARAM(PARAM_GR_I_IL10R_TCELL));
}

void Tcell::moveTcell(GrGrid& grid, bool ccl2, bool ccl5, bool cxcl9)
{
  Pos pos  = Agent::moveAgent(grid, ccl2, ccl5, cxcl9, false, _PARAM(PARAM_TCELL_MOVEMENT_BONUSFACTOR));

  // Check whether newCell is not caseated and contains empty slots
  size_t n = grid.getNumberOfAgents(pos);
  if (n != 2 && !grid.isCaseated(pos))
    {
      // Move with p = 1, if newCell is empty
      // Move with p = _PROB_MOVE_TCELL_TO_MAC, if newCell contains a macrophage
      // Move with p = _PROB_MOVE_TCELL_TO_TCELL, if newCell contains a T cell
      if ((n == 0) ||
          (grid.hasAgentType(MAC, pos) && g_Rand.getReal() < _PARAM(PARAM_TCELL_PROB_MOVE_TO_MAC)) ||
          (grid.hasTcell(pos) && g_Rand.getReal() < _PARAM(PARAM_TCELL_PROB_MOVE_TO_TCELL)))
        {
          assert_res(grid.removeAgent(this));
          grid.addAgent(this, pos);
          _pos = pos;
        }
    }
}
