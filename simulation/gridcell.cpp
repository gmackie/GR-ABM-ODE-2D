/*
 * gridcell.cpp
 *
 *  Created on: 02-nov-2009
 *      Author: M. El-Kebir
 */

#include "gridcell.h"
#include "grgrid.h"
#include "serialization.h"

const std::string GridCell::_ClassName = "GridCell";

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

	Serialization::writeHeader(out, GridCell::_ClassName);

	out << _row << std::endl;
	out << _col << std::endl;
	out << _source << std::endl;
	out << _nKillings << std::endl;
	out << _nRecruitments << std::endl;
	out << _nSecretions << std::endl;
	out << _macAttractant << std::endl;
	out << _tnf << std::endl;
	out << _ccl2 << std::endl;
	out << _ccl5 << std::endl;
	out << _cxcl9 << std::endl;
	out << _shedtnfr2 << std::endl;
	out << _extMtb << std::endl;

	Serialization::writeFooter(out, GridCell::_ClassName);
}

void GridCell::deserialize(std::istream& in)
{
	assert(in.good());

	if (!Serialization::readHeader(in, GridCell::_ClassName))
	{
		exit(1);
	}

	in >> _row;
	in >> _col;
	in >> _source;
	in >> _nKillings;
	in >> _nRecruitments;
	in >> _nSecretions;
	in >> _macAttractant;
	in >> _tnf;
	in >> _ccl2;
	in >> _ccl5;
	in >> _cxcl9;
	in >> _shedtnfr2;
	in >> _extMtb;

	if (!Serialization::readFooter(in, GridCell::_ClassName))
	{
		exit(1);
	}

	// These get assigned when the agent lists are deserialized and processed.
	_agent[0] = _agent[1] = NULL;
}

int GridCell::getOccupiedNeighborCount(const GrGrid& grid) const
{
	int mcCount = 0; // The number of compartments in the Moore neighborhood that are occupied.

	// Check the micro-compartment's Moore neighborhood (including the micro-compartment itself).
	for (int i = -1; i <= 1; i++)
	{
		for (int j = -1; j <= 1; j++)
		{
			if (grid(MOD_ROW(_row + i), MOD_COL(_col + j)).isOccupied())
			{
				mcCount++;
			}
		}
	}

	return mcCount;
}

float GridCell::getCellDensity(const GrGrid& grid) const
{
	float cellDensity = 0.0; // The cell density in the Moore neighborhood of this micro-compartment.
	float mcCount = getOccupiedNeighborCount(grid); // The number of compartments in the Moore neighborhood that are occupied.

	cellDensity = mcCount/MOORE_COUNT_DBL;

	return cellDensity;
}
