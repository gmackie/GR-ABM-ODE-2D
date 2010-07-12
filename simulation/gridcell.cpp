/*
 * gridcell.cpp
 *
 *  Created on: 02-nov-2009
 *      Author: M. El-Kebir
 */

#include "gridcell.h"

GridCell::GridCell()
	: _row(-1)
	, _col(-1)
	, _source(false)
	, _nKillings(0)
	, _nRecruitments(0)
	, _nSecretions(0)
	, _macAttractant(0)
	, _tnf(0)
	, _ccl2(0)
	, _ccl5(0)
	, _cxcl9(0)
	, _shedtnfr2(0)
	, _extMtb(0)
{
	_agent[0] = _agent[1] = NULL;
}

GridCell::~GridCell()
{
}

bool GridCell::addAgent(Agent* pAgent)
{
	// check if there is space left
	if (_agent[0] && _agent[1])
		return false;

	// check if cell is not caseated
	if (isCaseated())
		return false;

	// we now have that either _agent[0] or _agent[1] is NULL
	// in other words: there is an empty slot

	if (Mac::isMac(pAgent) && hasMac())
	{
		// at most one macrophage can be present in a grid cell
		return false;
	}
	else 
	{
		// put pAgent in an empty slot
		if (_agent[0])
		{
			_agent[1] = pAgent;
		}
		else
		{
			_agent[0] = pAgent;
		}
		return true;
	}
}

void GridCell::serialize(std::ostream& out) const
{
	assert(out.good());

	out << _row << std::endl;
	out << _col << std::endl;
	out << _source << std::endl;
	out << _nKillings << std::endl;
	out << _nRecruitments << std::endl;
	out << _tnf << std::endl;
	out << _ccl2 << std::endl;
	out << _ccl5 << std::endl;
	out << _cxcl9 << std::endl;
	out << _shedtnfr2 << std::endl;
	out << _extMtb << std::endl;
	out << _macAttractant << std::endl;
}

void GridCell::deserialize(std::istream& in)
{
	assert(in.good());

	in >> _row;
	in >> _col;
	in >> _source;
	in >> _nKillings;
	in >> _nRecruitments;
	in >> _tnf;
	in >> _ccl2;
	in >> _ccl5;
	in >> _cxcl9;
	in >> _shedtnfr2;
	in >> _extMtb;
	in >> _macAttractant;

	_agent[0] = _agent[1] = NULL;
}
