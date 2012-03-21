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

unsigned long Agent::nextID = 0;

// Needed for deserializing the model state.
// Avoids the calls to the random number generator in the normal constructor, allowing the random number generator
// to remain in synch after deserialization.
Agent::Agent()
	: _id(0)
	, _birthTime(-1)
	, _deathTime(-1)
	, _pos(-1, -1)

	// TNFR components
	, _mTNF(-1.0)
	, _surfTNFR1(-1.0)
	, _surfTNFR2(-1.0)
	, _surfBoundTNFR1(-1.0)
	, _surfBoundTNFR2(-1.0)
	, _intBoundTNFR1(-1.0)
	, _intBoundTNFR2(-1.0)
	, _mTNFRNA(-1.0)
	, _vTNFR1(-1.0)
	, _vTNFR2(-1.0)
	, _kSynth(-1.0)
	, _kTACE(-1.0)
	, _kmRNA(-1.0)

	// IL10 components
	, _surfIL10R(-1.0)
	, _vIL10R(-1.0)
	, _surfBoundIL10R(-1.0)
	, _kISynth(-1.0)

{
}

Agent::Agent(int birthtime, int deathtime, int row, int col

				//TNFR Components
				, Scalar meanTNFR1
				, Scalar stdTNFR1
				, Scalar meanTNFR2
				, Scalar stdTNFR2
				, Scalar kSynth
				, Scalar kTACE

				// IL10 components
				, Scalar iIL10R
				, Scalar stdIL10R
			)
	: _id(0)
	, _birthTime(birthtime)
	, _deathTime(deathtime)
	, _pos(row, col)

	// TNFR components
	, _mTNF(0.0)
	, _surfTNFR1(g_Rand.getReal(meanTNFR1 - stdTNFR1, meanTNFR1 + stdTNFR1))
	, _surfTNFR2(g_Rand.getReal(meanTNFR2 - stdTNFR2, meanTNFR2 + stdTNFR2))
	//	, _surfTNFR1(g_Rand.getLogNormal(meanTNFR1),_PARAM(stdTNFR1)))
	//	, _surfTNFR2(g_Rand.getLogNormal(meanTNFR2),_PARAM(stdTNFR2)))
	, _surfBoundTNFR1(0.0)
	, _surfBoundTNFR2(0.0)
	, _intBoundTNFR1(0.0)
	, _intBoundTNFR2(0.0)
	, _mTNFRNA(0.0)
	, _vTNFR1(_surfTNFR1 * _PARAM(PARAM_GR_K_T1))
	, _vTNFR2(_surfTNFR2 * _PARAM(PARAM_GR_K_T2))
	, _kSynth(kSynth)
	, _kTACE(kTACE)
	, _kmRNA(0.0)

	// IL10 components
	, _surfIL10R(g_Rand.getReal(iIL10R - stdIL10R, iIL10R + stdIL10R))
	, _vIL10R(_surfIL10R * _PARAM(PARAM_GR_I_K_T))
	, _surfBoundIL10R(0.0)
	, _kISynth(0.0)
{
	_id = createID();
}

Agent::~Agent()
{
}

Pos Agent::moveAgent(GrGrid& grid, bool ccl2, bool ccl5, bool cxcl9, bool attractant, double bonusFactor)
{
	int DestinationOrdinal = getDestinationOrdinal(grid, ccl2, ccl5, cxcl9, attractant, bonusFactor);
	return compartmentOrdinalToCoordinates(DestinationOrdinal, grid.getRange());
}

int Agent::getDestinationOrdinal(GrGrid& grid, bool ccl2, bool ccl5, bool cxcl9, bool attractant, double bonusFactor)
{

	double min = _PARAM(PARAM_GR_MIN_CHEMOTAXIS);
	double max = _PARAM(PARAM_GR_MAX_CHEMOTAXIS);

	bool ccl2Switch =  min < grid.CCL2(_pos) && grid.CCL2(_pos) < max;
	bool ccl5Switch =  min < grid.CCL5(_pos) && grid.CCL5(_pos) < max;
	bool cxcl9Switch = min < grid.CXCL9(_pos)  && grid.CXCL9(_pos) < max;

	Scalar prob[MOORE_COUNT] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

  int _row(GETROW(_pos)), _col(GETCOL(_pos));
	int k = 0;
  for (int i = -1; i <= 1; i++)
  {
    for (int j = -1; j <= 1; j++)
    {
      if (ccl2 && ccl2Switch)
        prob[k] += grid.CCL2(Pos(grid.mod_row(_row + i), grid.mod_col(_col + j)));
      if (ccl5 && ccl5Switch)
        prob[k] += grid.CCL5(Pos(grid.mod_row(_row + i), grid.mod_col(_col + j)));
      if (cxcl9 && cxcl9Switch)
        prob[k] += grid.CXCL9(Pos(grid.mod_row(_row + i), grid.mod_col(_col + j)));
      if (attractant)
        prob[k] += grid.macAttractant(Pos(grid.mod_row(_row + i), grid.mod_col(_col + j)));
      k++;
    }
  }
	
	// multiply highest probability with bonusFactor,
	// and while we are at it, determine the sum
	k = 0;
	double sum = 0;
	for (int i = 0; i < MOORE_COUNT; i++)
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
		for (int i = 0; i < MOORE_COUNT; i++)
			prob[i] /= sum;

		// compute cumulative array
		double cumProb[MOORE_COUNT];
    cumProb[0] = prob[0];
		for (int i = 1; i < MOORE_COUNT; i++)
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
  assert(0 <= k && k < MOORE_COUNT);
	return k;
}

Pos Agent::compartmentOrdinalToCoordinates(int ordinal, const Pos& dim) const
{

	/**
	 * The possible values for the ordinal are:
	 *
	 *  0  1  2
	 *  3  4  5
	 *  6  7  8
	 */
	int dRow = ((ordinal / 3) % 3) - 1;
	int dCol = ordinal % 3 - 1;
	int newRow = (_pos.x + dRow + dim.x) % dim.x;
	int newCol = (_pos.y + dCol + dim.y) % dim.y;

	return Pos(newRow, newCol);

}

void Agent::serialize(std::ostream& out) const
{
	assert(out.good());

	Serialization::writeHeader(out, Agent::_ClassName);

	out << _birthTime << std::endl;
	out << _deathTime << std::endl;
	out << _pos.x << std::endl;
	out << _pos.y << std::endl;

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
	in >> _pos.x;
	in >> _pos.y;

	if (!Serialization::readFooter(in, Agent::_ClassName))
	{
		exit(1);
	}
}
