/*
 * agent.cpp
 *
 *  Created on: 02-nov-2009
 *      Author: M. El-Kebir
 */

#include "agent.h"
#include "grgrid.h"

Agent::Agent(int birthtime, int deathtime, int row, int col)
	: _birthTime(birthtime)
	, _deathTime(deathtime)
	, _row(row)
	, _col(col)
{
	assert(0 <= _row && _row < NROWS);
	assert(0 <= _col && _col < NCOLS);
}

Agent::~Agent()
{
}

int Agent::moveAgent(GrGrid& grid, bool ccl2, bool ccl5, bool cxcl9, bool attractant, double bonusFactor)
{
	GridCell& cell = grid(_row, _col);

	double min = _PARAM(PARAM_GR_MIN_CHEMOTAXIS);
	double max = _PARAM(PARAM_GR_MAX_CHEMOTAXIS);

	bool ccl2Switch = min < cell.getCCL2() && cell.getCCL2() < max;
	bool ccl5Switch = min < cell.getCCL5() && cell.getCCL5() < max;
	bool cxcl9Switch = min < cell.getCXCL9() && cell.getCXCL9() < max;

	double prob[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

	int k = 0;
	for (int i = -1; i <= 1; i++)
	{
		for (int j = -1; j <= 1; j++)
		{
			if (ccl2 && ccl2Switch)
				prob[k] += grid(MOD_ROW(_row + i), MOD_COL(_col + j)).getCCL2();
			if (ccl5 && ccl5Switch)
				prob[k] += grid(MOD_ROW(_row + i), MOD_COL(_col + j)).getCCL5();
			if (cxcl9 && cxcl9Switch)
				prob[k] += grid(MOD_ROW(_row + i), MOD_COL(_col + j)).getCXCL9();
			if (attractant)
				prob[k] += grid(MOD_ROW(_row + i), MOD_COL(_col + j)).getMacAttractant();

			// add 1 just to prevent that cumArray contains only zeros
			prob[k++]++;
		}
	}
	
	// multiply highest probability with bonusFactor,
	// and while we are at it, determine the sum
	k = 0;
	double sum = 0;
	for (int i = 0; i < 9; i++)
	{
		if (prob[k] > prob[i])
			k = i;
		sum += prob[i];
	}
	sum += (bonusFactor - 1) * prob[k];
	prob[k] *= bonusFactor;

	// normalize
	for (int i = 0; i < 9; i++)
		prob[i] /= sum;

	// compute cumulative array
	double cumProb[9] = {prob[0], 0, 0, 0, 0, 0, 0, 0, 0};
	for (int i = 1; i < 9; i++)
	{
		cumProb[i] = (cumProb[i - 1] + prob[i]);
	}

	// linear search
	double r = g_Rand.getReal();
	for (k = 0; cumProb[k] < r && k < 9; k++);

	/**
	 * 0 1 2
	 * 3 4 5
	 * 6 7 8
	 */
	return k;
}

void Agent::serialize(std::ostream& out) const
{
	assert(out.good());

	out << _birthTime << std::endl;
	out << _deathTime << std::endl;
	out << _row << std::endl;
	out << _col << std::endl;
}

void Agent::deserialize(std::istream& in)
{
	assert(in.good());

	in >> _birthTime;
	in >> _deathTime;
	in >> _row;
	in >> _col;
}
