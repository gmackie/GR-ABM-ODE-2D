/*
 * grid.cpp
 *
 *  Created on: 02-nov-2009
 *      Author: M. El-Kebir
 */

#include "grgrid.h"
#include "serialization.h"
#include "agent.h"

using namespace std;

const std::string GrGrid::_ClassName = "GrGrid";


#define SZ (dim.x*dim.y)
#define Scalar_SZ ((dim.x+2)*(dim.y+2))
#define GRID(type, name) \
  , _##name (SZ)
#define PADDED_GRID(type, name) \
  , _##name(Scalar_SZ)	\
  , _u_##name(Scalar_SZ)	\
  , _v_##name(Scalar_SZ)

GrGrid::GrGrid(const Pos& dim)
  : _dim(dim)
  , _nCaseation(_PARAM(PARAM_GR_NR_KILLINGS_FOR_CASEATION))
  , _sources()
  , _agents(SZ*MAX_AGENTS_PER_CELL, NULL)
  GRIDS_DEFS
{
  
}
#undef GRID
#undef PADDED_GRID
#undef SZ
#undef Scalar_SZ

GrGrid::~GrGrid()
{
}

void GrGrid::initSources()
{
	int nSources = _PARAM(PARAM_GR_NR_SOURCES);

	// Useful for testing and debugging.
	// Prevents recruitment and can just watch initial macs without the sources cluttering the screen.
	if (nSources == 0)
		return;
  if (nSources > getSize())
  {
    std::cerr <<"Cannot place "<<nSources<<" on grid "<<_dim.x<<'x'<<_dim.y<<std::endl;
    exit(1);
  }

	const int n = (int)floor(sqrt((float) nSources));

	const int dRow = _dim.x / n;
	const int dCol = _dim.y / n;
  const int leftover = _dim.x - n*dRow;  //Assumes square grid

  bool sources[_dim.x][_dim.y];
  memset(sources, 0, getSize()*sizeof(bool));

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			const int row = g_Rand.getInt(dRow);
			const int col = g_Rand.getInt(dCol);
      const Pos p(dRow*i+row+std::min(i, leftover), dCol*j+col+std::min(j, leftover));

			assert(!sources[p.x][p.y]);
      sources[p.x][p.y] = true;
			_sources.push_back(p);
		}
	}

	// pick remaining sources
	while ((int)_sources.size() < nSources)
	{
		const int row = g_Rand.getInt(_dim.x);
		const int col = g_Rand.getInt(_dim.y);

		if (!sources[row][col])
		{
      sources[row][col] = true;
			_sources.push_back(Pos(row, col));
		}
	}
}

void GrGrid::shuffleSources()
{
	random_shuffle(_sources.begin(), _sources.end(), g_Rand);
}

void GrGrid::serialize(std::ostream& out) const
{
	assert(out.good());

	Serialization::writeHeader(out, GrGrid::_ClassName);

	out << getRange() << endl;

	out << _nCaseation << endl;

	Pos p;

	// The main grid.
	for (p.x = 0; p.x < _dim.x; p.x++)
	{
		for (p.y = 0; p.y < _dim.y; p.y++)
		{
			out << nKillings(p) << endl;
            out << trappedCaseation(p) << endl;
			out << nRecruitments(p) << endl;
            out << nRecruitmentsMac(p) << endl;
            out << nRecruitmentsTgam(p) << endl;
            out << nRecruitmentsTcyt(p) << endl;
            out << nRecruitmentsTreg(p) << endl;
			out << nSecretions(p) << endl;
			out << macAttractant(p) << endl;
			out << TNF(p) << endl;
			out << CCL2(p) << endl; // CCL5, CSCL9 values are based on CCL2 values.
			out << shedTNFR2(p) << endl;
			out << il10(p) << endl;
			out << extMTB(p) << endl;
			out << std::endl;
		}
		out << std::endl;
	}

	// The u grid.
	for (p.x = 0; p.x < _dim.x; p.x++)
	{
		for (p.y = 0; p.y < _dim.y; p.y++)
		{
			out << u_macAttractant(p) << endl;
			out << u_TNF(p) << endl;
			out << u_CCL2(p) << endl; // CCL5, CSCL9 values are based on CCL2 values.
			out << u_shedTNFR2(p) << endl;
			out << u_il10(p) << endl;
			out << std::endl;
		}
		out << std::endl;
	}

	// The v grid.
	for (p.x = 0; p.x < _dim.x; p.x++)
	{
		for (p.y = 0; p.y < _dim.y; p.y++)
		{
			out << v_macAttractant(p) << endl;
			out << v_TNF(p) << endl;
			out << v_CCL2(p) << endl; // CCL5, CSCL9 values are based on CCL2 values.
			out << v_shedTNFR2(p) << endl;
			out << v_il10(p) << endl;
			out << std::endl;
		}
		out << std::endl;
	}

	out <<_sources.size()<<std::endl;
	for(std::vector<Pos>::const_iterator v= _sources.begin(); v!=_sources.end(); v++)
	{
		out << *v << endl;
	}

	Serialization::writeFooter(out, GrGrid::_ClassName);
}

void GrGrid::deserialize(std::istream& in)
{
  assert(in.good());

  Serialization::readHeader(in, GrGrid::_ClassName);

  Pos range;
  const Pos& standard = getRange();
  in>>range;

  if(GETROW(range) != GETROW(standard) && GETCOL(range) != GETCOL(standard))
    throw std::length_error("Dimension mismatch in deserialization");

  in >> _nCaseation;

  _sources.clear();

  Pos p;

  // The main grid.
  for(p.x = 0; p.x < _dim.x; p.x++)
    for(p.y = 0; p.y < _dim.y; p.y++)
    {

		in >> (int&) nKillings(p);
        in >> (int&) trappedCaseation(p);
		in >> (int&) nRecruitments(p);
        in >> (int&) nRecruitmentsMac(p);
        in >> (int&) nRecruitmentsTgam(p);
        in >> (int&) nRecruitmentsTcyt(p);
        in >> (int&) nRecruitmentsTreg(p);
		in >> (int&) nSecretions(p);
		in >> (Scalar&) macAttractant(p);
		in >> (Scalar&) TNF(p);
		in >> (Scalar&) CCL2(p);  // CCL5, CSCL9 values are based on CCL2 values.
		in >> (Scalar&) shedTNFR2(p);
		in >> (Scalar&) il10(p);
		in >> (Scalar&) extMTB(p);

        // Clear the agents, in case some had been defined during simulation initialization,
        // for example initial agents defined from a parameter file INIT section.
        for (size_t i = 0; i < MAX_AGENTS_PER_CELL; i++)
        {
        	agent(p, i) = NULL;
        }
     }

  // The u grid.
  for(p.x = 0; p.x < _dim.x; p.x++)
    for(p.y = 0; p.y < _dim.y; p.y++)
    {
		in >> (Scalar&) u_macAttractant(p);
		in >> (Scalar&) u_TNF(p);
		in >> (Scalar&) u_CCL2(p);  // CCL5, CSCL9 values are based on CCL2 values.
		in >> (Scalar&) u_shedTNFR2(p);
		in >> (Scalar&) u_il10(p);
     }

  // The v grid.
  for(p.x = 0; p.x < _dim.x; p.x++)
    for(p.y = 0; p.y < _dim.y; p.y++)
    {
		in >> (Scalar&) v_macAttractant(p);
		in >> (Scalar&) v_TNF(p);
		in >> (Scalar&) v_CCL2(p);  // CCL5, CSCL9 values are based on CCL2 values.
		in >> (Scalar&) v_shedTNFR2(p);
		in >> (Scalar&) v_il10(p);
     }

  int sz = 0;
  in >> sz;
  for(int i=0;i<sz;i++) {
    in>>p;
    _sources.push_back(p);
  }

  Serialization::readFooter(in, GrGrid::_ClassName);
}

bool GrGrid::incKillings(const Pos& p) {
  assert(inRange(p));
  ++_nKillings[Indexer::ind(_dim, p)];
  if(isCaseated(p)) {
    for(unsigned i=0;i<MAX_AGENTS_PER_CELL;i++)
      if(agent(p, i) && !agent(p,i)->isDead())
      {
          agent(p,i)->kill();
          ++_nKillings[Indexer::ind(_dim, p)];
      }
    return true;
  }
  return false;
}

int GrGrid::isTrapped(const Pos& p)
{
    assert(inRange(p));
    int caseationCount = 0;
    for (int i=-1; i<=1; i++)
    {
        for (int j=-1; j<=1; j++)
        {
            Pos modp(mod_row(p.x + i), mod_col(p.y + j));
            if (isCaseated(modp))
                caseationCount++;
        }
    }

    if (caseationCount == (MOORE_COUNT - 1))
    {
        ++_trappedCaseation[Indexer::ind(_dim,p)];
//        std::cout << "Trapped Cell to Caseation @: " << p << std::endl;
        // Kill any cell that resides in the trapped compartment
        for(unsigned i=0;i<MAX_AGENTS_PER_CELL;i++)
          if(agent(p, i) && !agent(p,i)->isDead())
          {
              agent(p,i)->kill();
              ++_nKillings[Indexer::ind(_dim, p)];
          }
    }
    return caseationCount;
}

int GrGrid::getOccupiedNeighborCount(int _row, int _col) const 
{ 
  int mcCount = 0; // The number of compartments in the Moore neighborhood that are occupied. 
 
  // Check the micro-compartment's Moore neighborhood (including the micro-compartment itself). 
  for (int k = -1; k <= 1; k++) 
  { 
    for (int i = -1; i <= 1; i++) 
    { 
      if (isOccupied(mod_row(_row + k), mod_col(_col + i))) 
      { 
        mcCount++; 
      } 
    } 
  } 
 
  return mcCount; 
} 

int GrGrid::getOccupiedNeighborCount(const Pos& p) const {
  return getOccupiedNeighborCount(GETROW(p), GETCOL(p));
}

Scalar GrGrid::getCellDensity(const Pos& p) const
{
  return getCellDensity(p.x, p.y);
}

Scalar GrGrid::getCellDensity(int _row, int _col) const
{
  Scalar cellDensity = 0.0; // The cell density in the Moore neighborhood of this micro-compartment.
  Scalar mcCount = getOccupiedNeighborCount(_row, _col); // The number of compartments in the Moore neighborhood that are occupied.

  cellDensity = mcCount/MOORE_COUNT_DBL;

  return cellDensity;
}

int GrGrid::hasAgentType(AgentType TYPE, const Pos& p) const {
  int sum = 0;
  for(unsigned i=0;i<MAX_AGENTS_PER_CELL;i++)
    if(agent(p,i) && agent(p,i)->getAgentType() == TYPE) sum++;
  return sum;
}

int GrGrid::hasTcell(const Pos& p) const {
  int sum = 0;
  for(unsigned i=0;i<MAX_AGENTS_PER_CELL;i++)
    if(agent(p,i) && agent(p,i)->getAgentType() != MAC) sum++;
  return sum;
}

int GrGrid::hasAgentType(AgentType TYPE, int s, const Pos& p) const {
  int sum = 0;
  for(unsigned i=0;i<MAX_AGENTS_PER_CELL;i++)
  {
    if(agent(p,i) && agent(p,i)->getAgentType() == TYPE && agent(p,i)->getState() == s) sum++;
  }
  return sum;
}

bool GrGrid::addAgent(Agent* a, const Pos& p) {
  assert(a);
  if(isCaseated(p)) return false;
  if(a->getAgentType() == MAC && hasAgentType(MAC, p)) return false;
  for(unsigned i=0;i<MAX_AGENTS_PER_CELL;i++){
    if(agent(p,i)) continue;
    agent(p,i) = a;
    return true;
  }
  return false;
}

bool GrGrid::addAgent(Agent* a, int i, int j) {
  return addAgent(a, Pos(i,j));
}

bool GrGrid::addAgent(Agent* a) {
  return addAgent(a, a->getPosition());
}

bool GrGrid::addAgent(Agent* a, const Pos& p, int index)
{
  assert(a);
  if(isCaseated(p)) return false;
  if(a->getAgentType() == MAC && hasAgentType(MAC, p)) return false;

  if(!agent(p, index))
  {
    agent(p, index) = a;
    return true;
  }
  return false;
}

bool GrGrid::addAgent(Agent* a, int index)
{
  return addAgent(a, a->getPosition(), index);
}

bool GrGrid::addAgent(Agent* a, int i, int j, int index)
{
  return addAgent(a, Pos(i, j), index);
}

bool GrGrid::removeAgent(const Agent* a) {
  const Pos& p = a->getPosition();
  for(unsigned i=0;i<MAX_AGENTS_PER_CELL;i++)
    if(a == agent(p,i)) {
      agent(p,i) = NULL;
      return true;
    }
  return false;
}

int GrGrid::agentIndex (const Agent* pAgent) const
{
        assert(pAgent);

        for(size_t i = 0 ; i < MAX_AGENTS_PER_CELL; i++)
        {
                const Agent* pGridAgent = agent(pAgent->getPosition(), i);
                if (pGridAgent == pAgent)
                {
                        return i;
                }
        }

        return -1;
}
