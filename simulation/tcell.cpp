/*
 * tcell.cpp
 *
 *  Created on: 13-nov-2009
 *      Author: M. El-Kebir
 */


#include "tcell.h"
#include "grgrid.h"

Tcell::Tcell(int birthtime, int row, int col)
	: Agent(birthtime, birthtime + _PARAM(PARAM_TCELL_AGE), row, col)
{
}

Tcell::~Tcell()
{
}

void Tcell::moveTcell(GrGrid& grid, bool ccl2, bool ccl5, bool cxcl9)
{
	int k = Agent::moveAgent(grid, ccl2, ccl5, cxcl9, false, 1);
	
	/**
	 * 0 1 2
	 * 3 4 5
	 * 6 7 8
	 */
	int dRow = k / 3 - 1;
	int dCol = k % 3 - 1;
	int newRow = MOD_ROW(_row + dRow);
	int newCol = MOD_COL(_col + dCol);

	GridCell& cell = grid(_row, _col);
	GridCell& newCell = grid(newRow, newCol);

	// Check whether newCell is not caseated and contains empty slots
	if (newCell.getNumberOfAgents() != 2 && !newCell.isCaseated())
	{
		// Move with p = 1, if newCell is empty
		// Move with p = _PROB_MOVE_TCELL_TO_MAC, if newCell contains a macrophage
		// Move with p = _PROB_MOVE_TCELL_TO_TCELL, if newCell contains a T cell
		if ((newCell.getNumberOfAgents() == 0) ||
			(newCell.hasMac() && g_Rand.getReal() < _PARAM(PARAM_TCELL_PROB_MOVE_TO_MAC)) ||
			(newCell.hasTcell() && g_Rand.getReal() < _PARAM(PARAM_TCELL_PROB_MOVE_TO_TCELL)))
		{
			_row = newRow;
			_col = newCol;
			assert_res(cell.removeAgent(this));
			newCell.addAgent(this);
		}
	}
}

void Tcell::serialize(std::ostream& out) const
{
	Agent::serialize(out);
}

void Tcell::deserialize(std::istream& in)
{
	Agent::deserialize(in);
}
