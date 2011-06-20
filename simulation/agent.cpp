/*
 * agent.cpp
 *
 *  Created on: 02-nov-2009
 *      Author: M. El-Kebir
 */

#include "agent.h"
#include "grgrid.h"
#include "serialization.h"

const std::string Agent::_ClassName = "Agent";

// Needed for deserializing the model state.
// Avoids the calls to the random number generator in the normal constructor, allowing the random number generator
// to remain in synch after deserialization.
Agent::Agent()
	: _birthTime(-1)
	, _deathTime(-1)
	, _row(-1)
	, _col(-1)
{
}

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

			k++;
		}
	}
	
	// multiply highest probability with bonusFactor,
	// and while we are at it, determine the sum
	k = 0;
	double sum = 0;
	for (int i = 0; i < 9; i++)
	{
		if (prob[i] > prob[k])
			k = i;
		sum += prob[i];
	}
	sum += (bonusFactor - 1) * prob[k];
	prob[k] *= bonusFactor;

	if (sum > 0.0)
	{
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
		for (k = 0; k < 9 && cumProb[k] < r; k++);
	}
	else
	{
		// prob[i] = 0 for all i.
		// Pick from the neighbors with equal probability.
		k = g_Rand.getInt(MOORE_COUNT, 0);
	}

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

	Serialization::writeHeader(out, Agent::_ClassName);

	out << _birthTime << std::endl;
	out << _deathTime << std::endl;
	out << _row << std::endl;
	out << _col << std::endl;

	Serialization::writeFooter(out, Agent::_ClassName);
}

void Agent::deserialize(std::istream& in)
{
	assert(in.good());

	if (!Serialization::readHeader(in, Agent::_ClassName))
	{
		exit(1);
	}

	in >> _birthTime;
	in >> _deathTime;
	in >> _row;
	in >> _col;

	if (!Serialization::readFooter(in, Agent::_ClassName))
	{
		exit(1);
	}
}
