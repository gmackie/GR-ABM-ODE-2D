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

unsigned long Agent::_nextID = 0;

// Needed for deserializing the model state.
// Avoids the calls to the random number generator in the normal constructor, allowing the random number generator
// to remain in synch after deserialization.
Agent::Agent()
	: _id(0)
	, _birthTime(-1)
	, _deathTime(-1)
	, _pos(-1, -1)
	, _trackMolecularDynamics(false)

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
	, _trackMolecularDynamics(false)

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

void Agent::solveTNF(GrGrid& grid, double dt)
{
	double koff1 = _PARAM(PARAM_GR_K_ON1) * _PARAM(PARAM_GR_KD1);
	double koff2 = _PARAM(PARAM_GR_K_ON2) * _PARAM(PARAM_GR_KD2);

	double tnf = grid.TNF(_pos) / (NAV * VOL);
	double shedtnfr2 = grid.shedTNFR2(_pos) / (NAV * VOL);
    double il10 = grid.il10(_pos) /(NAV * VOL);

    double dmTNFRNA;
	double dmTNF;
	double dsurfTNFR1;
	double dsurfTNFR2;
	double dsurfBoundTNFR1;
	double dsurfBoundTNFR2;
	double dintBoundTNFR1;
	double dintBoundTNFR2;
	double dsTNF;
	double dshedTNFR2;

    double eqsurfBoundIL10R;
    double IkmRNA;

    // modulate TNF parameters based on equilibrium IL10 receptor equations since the il10 molecular equations are not active
    eqsurfBoundIL10R = (il10 * _surfIL10R) / (_PARAM(PARAM_GR_I_KD) + il10);

    if (_kmRNA > 0) {
        IkmRNA = _kmRNA * ((_kSynth/_kmRNA)+ ((1.0 - (_kSynth/_kmRNA))/(1.0 + pow(2.7183, ((eqsurfBoundIL10R - _PARAM(PARAM_GR_LINK_RNA_GAMMA))/_PARAM(PARAM_GR_LINK_RNA_DELTA))))));
    }
    else
    {
        IkmRNA = 0.0;
    }

    // end of equilibrium calculations

    dmTNFRNA = (IkmRNA - _PARAM(PARAM_GR_K_TRANS) * _mTNFRNA) * dt;
	dmTNF = (_PARAM(PARAM_GR_K_TRANS) * _mTNFRNA - _kTACE * _mTNF) * dt;
	dsurfTNFR1 = (_vTNFR1 - _PARAM(PARAM_GR_K_ON1) * tnf * _surfTNFR1 + koff1 * _surfBoundTNFR1 - _PARAM(PARAM_GR_K_T1) * _surfTNFR1 + _PARAM(PARAM_GR_K_REC1) * _intBoundTNFR1) * dt;
	dsurfTNFR2 = (_vTNFR2 - _PARAM(PARAM_GR_K_ON2) * tnf * _surfTNFR2 + koff2 * _surfBoundTNFR2 - _PARAM(PARAM_GR_K_T2) * _surfTNFR2 + _PARAM(PARAM_GR_K_REC2) * _intBoundTNFR2) * dt;
	dsurfBoundTNFR1 = (_PARAM(PARAM_GR_K_ON1) * tnf * _surfTNFR1 - koff1 * _surfBoundTNFR1 - _PARAM(PARAM_GR_K_INT1) * _surfBoundTNFR1) * dt;
	dsurfBoundTNFR2 = (_PARAM(PARAM_GR_K_ON2) * tnf * _surfTNFR2 - koff2 * _surfBoundTNFR2 - _PARAM(PARAM_GR_K_INT2) * _surfBoundTNFR2 - _PARAM(PARAM_GR_K_SHED) * _surfBoundTNFR2) * dt;
	dintBoundTNFR1 = (_PARAM(PARAM_GR_K_INT1) * _surfBoundTNFR1 - _PARAM(PARAM_GR_K_DEG1) * _intBoundTNFR1 - _PARAM(PARAM_GR_K_REC1) * _intBoundTNFR1) * dt;
	dintBoundTNFR2 = (_PARAM(PARAM_GR_K_INT2) * _surfBoundTNFR2 - _PARAM(PARAM_GR_K_DEG2) * _intBoundTNFR2 - _PARAM(PARAM_GR_K_REC2) * _intBoundTNFR2) * dt;
	dsTNF = ((DENSITY/NAV) * (_kTACE * _mTNF - _PARAM(PARAM_GR_K_ON1) * tnf * _surfTNFR1 + koff1 * _surfBoundTNFR1 - _PARAM(PARAM_GR_K_ON2) * tnf * _surfTNFR2 + koff2 * _surfBoundTNFR2)) * dt;
	dshedTNFR2 = ((DENSITY/NAV) * _PARAM(PARAM_GR_K_SHED) * _surfBoundTNFR2) * dt;

    _mTNFRNA += dmTNFRNA;
	_mTNF += dmTNF;
	_surfTNFR1 += dsurfTNFR1;
	_surfTNFR2 += dsurfTNFR2;
	_surfBoundTNFR1 += dsurfBoundTNFR1;
	_surfBoundTNFR2 += dsurfBoundTNFR2;
	_intBoundTNFR1 += dintBoundTNFR1;
	_intBoundTNFR2 += dintBoundTNFR2;
	tnf += dsTNF;
	shedtnfr2 += dshedTNFR2;

	grid.setTNF(_pos, (NAV * VOL * tnf));
	grid.setshedTNFR2(_pos, (NAV * VOL * shedtnfr2));
	if (_mTNF < 0 || _surfTNFR1 < 0 || _surfBoundTNFR1 < 0 || _surfTNFR2 < 0 || _surfBoundTNFR2 < 0 || _mTNFRNA < 0)
		std::cout << "Error: Negative Value of Species in TNF/TNFR dynamics" << std::endl;

    //cout << "Debug: Running TNF dynamics" << std::endl;
}

void Agent::solveIL10(GrGrid& grid, double dt)
{
    double il10 = grid.il10(_pos) / (NAV * VOL);

    double dsIL10;
	double dsurfIL10R;
	double dsurfBoundIL10R;

    // IL10 differential equations
    dsIL10 = (DENSITY/NAV) * _kISynth + ((DENSITY/NAV) * (_PARAM(PARAM_GR_I_K_OFF) * _surfBoundIL10R - _PARAM(PARAM_GR_I_K_ON) * _surfIL10R * il10)) * dt;
	dsurfIL10R = (_vIL10R - _PARAM(PARAM_GR_I_K_ON) * _surfIL10R * il10 + _PARAM(PARAM_GR_I_K_OFF) * _surfBoundIL10R - _PARAM(PARAM_GR_I_K_T) * _surfIL10R) * dt;
	dsurfBoundIL10R = (_PARAM(PARAM_GR_I_K_ON) * _surfIL10R * il10 - _PARAM(PARAM_GR_I_K_OFF) * _surfBoundIL10R - _PARAM(PARAM_GR_I_K_INT) * _surfBoundIL10R) * dt;
    // end of IL10 differential equations

    // update il10 variables
    _surfIL10R += dsurfIL10R;
    _surfBoundIL10R += dsurfBoundIL10R;
    il10 += dsIL10;

    grid.setil10(_pos, (NAV * VOL * il10));

    if (_surfIL10R < 0 || _surfBoundIL10R < 0)
        std::cout << "Error: Negative value of species in IL10/IL10R dynamics" << std::endl;

     //cout << "Debug: Running IL10 dynamics" << std::endl;

}

void Agent::solveDegradation(GrGrid& grid, double dt, bool tnfrDynamics, bool il10rDynamics, Scalar meanTNFR1, Scalar iIL10R)
{
    if (!tnfrDynamics) {

        // simulate the effect of TNF internalization by cells in the form of degradation. Only for TNF
        double dtnf;
        double tnf = grid.TNF(_pos);
        dtnf = -_PARAM(PARAM_GR_K_INT1) * (tnf / (tnf + _PARAM(PARAM_GR_KD1) * NAV * VOL)) * meanTNFR1 * dt * 0.4;
        grid.incTNF(_pos, dtnf);
    }

    if (!il10rDynamics) {

        double dil10;
        double il10 = grid.il10(_pos);

        // simulate the effect of IL10 internalization in the form of degradation. Only for IL10
        dil10 = -_PARAM(PARAM_GR_I_K_INT) * (il10 / (il10 + _PARAM(PARAM_GR_I_KD) * NAV * VOL)) * iIL10R * dt * _PARAM(PARAM_GR_I_MOD);
        grid.incil10(_pos, dil10);

    }

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


void Agent::classSerialize(std::ostream& out)
{
	assert(out.good());
	Serialization::writeHeader(out, Agent::_ClassName);
	out << _nextID << std::endl;
	Serialization::writeFooter(out, Agent::_ClassName);
}

void Agent::classDeserialize(std::istream& in)
{
	assert(in.good());

	if (!Serialization::readHeader(in, Agent::_ClassName))
	{
		exit(1);
	}

	in >> _nextID;

	if (!Serialization::readFooter(in, Agent::_ClassName))
	{
		exit(1);
	}
}


void Agent::serialize(std::ostream& out) const
{
	assert(out.good());

	Serialization::writeHeader(out, Agent::_ClassName);

	out << _id << std::endl;
	out << _birthTime << std::endl;
	out << _deathTime << std::endl;
	out << _pos.x << std::endl;
	out << _pos.y << std::endl;
	out << _trackMolecularDynamics << std::endl;

	// TNF associated attributes
	out << _mTNF << std::endl;
	out << _surfTNFR1 << std::endl;
	out << _surfTNFR2 << std::endl;
	out << _surfBoundTNFR1 << std::endl;
	out << _surfBoundTNFR2 << std::endl;
	out << _intBoundTNFR1 << std::endl;
	out << _intBoundTNFR2 << std::endl;
    out << _mTNFRNA << std::endl;
	out << _vTNFR1 << std::endl;
	out << _vTNFR2 << std::endl;
	out << _kSynth << std::endl;
	out << _kTACE << std::endl;
    out << _kmRNA << std::endl;

    // IL10 associated attributes
    out << _surfIL10R << std::endl;
    out << _vIL10R << std::endl;
    out << _surfBoundIL10R << std::endl;
    out << _kISynth << std::endl;

	Serialization::writeFooter(out, Agent::_ClassName);
}

void Agent::deserialize(std::istream& in)
{
	assert(in.good());

	if (!Serialization::readHeader(in, Agent::_ClassName))
	{
		exit(1);
	}

	in >> _id;
	in >> _birthTime;
	in >> _deathTime;
	in >> _pos.x;
	in >> _pos.y;
	in >> _trackMolecularDynamics;

	// TNF associated attributes
	in >> _mTNF;
	in >> _surfTNFR1;
	in >> _surfTNFR2;
	in >> _surfBoundTNFR1;
	in >> _surfBoundTNFR2;
	in >> _intBoundTNFR1;
	in >> _intBoundTNFR2;
    in >> _mTNFRNA;
	in >> _vTNFR1;
	in >> _vTNFR2;
	in >> _kSynth;
	in >> _kTACE;
    in >> _kmRNA;

    // IL10 associated attributes
    in >> _surfIL10R;
    in >> _vIL10R;
    in >> _surfBoundIL10R;
    in >> _kISynth;

	if (!Serialization::readFooter(in, Agent::_ClassName))
	{
		exit(1);
	}
}
