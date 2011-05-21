/*
 * gridcell.h
 *
 *  Created on: 02-nov-2009
 *      Author: M. El-Kebir
 */

#ifndef GRIDCELL_H
#define GRIDCELL_H

#include "gr.h"
#include "params.h"
#include "agent.h"
#include "tcell.h"
#include "macrophage.h"
#include "tregulatory.h"
#include "tcytotoxic.h"
#include "tgamma.h"

class GridCell
{
private:
	int _row;
	int _col;
	bool _source;
	int _nKillings;
	int _nRecruitments;
	int _nSecretions;
	double _macAttractant;
	double _tnf;
	double _ccl2;
	double _ccl5;
	double _cxcl9;
	double _shedtnfr2;
	double _extMtb;
	Agent* _agent[2];

public:
	GridCell();
	~GridCell();
	int getRow() const;
	int getCol() const;
	void setRowCol(int row, int col);
	double getMacAttractant() const;
	void setMacAttractant(double macAttractant);
	void incMacAttractant(double dMacAttractant);
	double getTNF() const;
	void setTNF(double tnf);
	void incTNF(double dTNF);
	double getCCL2() const;
	void setCCL2(double ccl2);
	void incCCL2(double dCCL2);
	double getCCL5() const;
	void setCCL5(double ccl5);
	void incCCL5(double dCCL5);
	double getCXCL9() const;
	void setCXCL9(double cxcl9);
	void incCXCL9(double dCXCL9);
	double getShedTNFR2() const;
	void setShedTNFR2(double shedtnfr2);
	void incShedTNFR2(double dShedTNFR2);
	double getExtMtb() const;
	void setExtMtb(double extMtb);
	void incExtMtb(double dExtMtb);
	bool isSource() const;
	void setSource();
	int getNrKillings() const;
	bool incNrKillings();
	int getNrRecruitments() const;
	void incNrRecruitments();
	int getNrSecretions() const;
	void incNrSecretions();
	bool isCaseated() const;
	bool addAgent(Agent* pAgent);
	bool hasMac() const;
	int hasTreg() const;
	int hasTcyt() const;
	int hasTgam() const;
	bool hasMac(MacState state) const;
	int hasTreg(TregState state) const;
	int hasTcyt(TcytState state) const;
	int hasTgam(TgamState state) const;
	int hasTcell() const;
	Agent* getAgent(int i);
	const Agent* getAgent(int i) const;
	bool removeAgent(Agent* pAgent);
	int getNumberOfAgents() const;
	void serialize(std::ostream& out) const;
	void deserialize(std::istream& in);
	bool isOccupied() const;
};

inline int GridCell::getRow() const
{
	return _row;
}

inline int GridCell::getCol() const
{
	return _col;
}

inline void GridCell::setRowCol(int row, int col)
{
	assert(0 <= row && row < NROWS);
	assert(0 <= col && col < NCOLS);

	_row = row;
	_col = col;
}

inline double GridCell::getMacAttractant() const
{
	return _macAttractant;
}

inline void GridCell::setMacAttractant(double macAttractant)
{
	_macAttractant = macAttractant;
}

inline void GridCell::incMacAttractant(double dMacAttractant)
{
	_macAttractant += dMacAttractant;
}

inline double GridCell::getTNF() const
{
	return _tnf;
}

inline void GridCell::setTNF(double tnf)
{
	_tnf = tnf;
}

inline void GridCell::incTNF(double dTNF)
{
	_tnf += dTNF;
}

inline double GridCell::getCCL2() const
{
	return _ccl2;
}

inline void GridCell::setCCL2(double ccl2)
{
	_ccl2 = ccl2;
}

inline void GridCell::incCCL2(double dCCL2)
{
	_ccl2 += dCCL2;
}

inline double GridCell::getCCL5() const
{
	return _ccl5;
}

inline void GridCell::setCCL5(double ccl5)
{
	_ccl5 = ccl5;
}

inline void GridCell::incCCL5(double dCCL5)
{
	_ccl5 += dCCL5;
}

inline double GridCell::getCXCL9() const
{
	return _cxcl9;
}

inline void GridCell::setCXCL9(double cxcl9)
{
	_cxcl9 = cxcl9;
}

inline void GridCell::incCXCL9(double dCXCL9)
{
	_cxcl9 += dCXCL9;
}

inline double GridCell::getShedTNFR2() const
{
	return _shedtnfr2;
}

inline void GridCell::setShedTNFR2(double shedtnfr2)
{
	_shedtnfr2 = shedtnfr2;
}

inline void GridCell::incShedTNFR2(double dShedTNFR2)
{
	_shedtnfr2 += dShedTNFR2;
}

inline double GridCell::getExtMtb() const
{
	return isCaseated() ? 0 : _extMtb;
}

inline void GridCell::setExtMtb(double extMtb)
{
	_extMtb = extMtb;
}

inline void GridCell::incExtMtb(double dExtMtb)
{
	_extMtb += dExtMtb;
}

inline bool GridCell::isSource() const
{
	return _source;
}

inline void GridCell::setSource()
{
	_source = true;
}

inline int GridCell::getNrKillings() const
{
	return _nKillings;
}

inline bool GridCell::incNrKillings()
{
	_nKillings++;
	
	if (isCaseated())
	{
		// clear T cells
		if (Tcell::isTcell(_agent[0]))
			_agent[0]->kill();

		if (Tcell::isTcell(_agent[1]))
			_agent[1]->kill();

		return true;
	}
	return false;
}

inline int GridCell::getNrRecruitments() const
{
	return _nRecruitments;
}

inline void GridCell::incNrRecruitments()
{
	_nRecruitments++;
}

inline int GridCell::getNrSecretions() const
{
	return _nSecretions;
}

inline void GridCell::incNrSecretions()
{
	_nSecretions++;
}

inline bool GridCell::isCaseated() const
{
	return _nKillings >= _PARAM(PARAM_GR_NR_KILLINGS_FOR_CASEATION);
}

inline bool GridCell::hasMac() const
{
	return Mac::isMac(_agent[0]) || Mac::isMac(_agent[1]);
}

inline int GridCell::hasTcell() const
{
	return (Tcell::isTcell(_agent[0]) ? 1 : 0) + (Tcell::isTcell(_agent[1]) ? 1 : 0);
}

inline int GridCell::hasTcyt() const
{
	return (Tcyt::isTcyt(_agent[0]) ? 1 : 0) + (Tcyt::isTcyt(_agent[1]) ? 1 : 0);
}

inline int GridCell::hasTgam() const
{
	return (Tgam::isTgam(_agent[0]) ? 1 : 0) + (Tgam::isTgam(_agent[1]) ? 1 : 0);
}

inline int GridCell::hasTreg() const
{
	return (Treg::isTreg(_agent[0]) ? 1 : 0) + (Treg::isTreg(_agent[1]) ? 1 : 0);
}

inline bool GridCell::hasMac(MacState state) const
{
	return Mac::isMac(_agent[0], state) || Mac::isMac(_agent[1], state);
}

inline int GridCell::hasTcyt(TcytState state) const
{
	return (Tcyt::isTcyt(_agent[0], state) ? 1 : 0) +
		(Tcyt::isTcyt(_agent[1], state) ? 1 : 0);
}

inline int GridCell::hasTgam(TgamState state) const
{
	return (Tgam::isTgam(_agent[0], state) ? 1 : 0) +
		(Tgam::isTgam(_agent[1], state) ? 1 : 0);
}

inline int GridCell::hasTreg(TregState state) const
{
	return (Treg::isTreg(_agent[0], state) ? 1 : 0) +
		(Treg::isTreg(_agent[1], state) ? 1 : 0);
}

inline Agent* GridCell::getAgent(int i)
{
	assert(0 <= i && i < 2);
	return _agent[i];
}

inline const Agent* GridCell::getAgent(int i) const
{
	assert(0 <= i && i < 2);
	return _agent[i];
}

inline bool GridCell::removeAgent(Agent* pAgent)
{
	for (int i = 0; i < 2; i++)
	{
		if (_agent[i] == pAgent)
		{
			_agent[i] = NULL;
			return true;
		}
	}
	return false;
}

inline int GridCell::getNumberOfAgents() const
{
	return (_agent[0] ? 1 : 0) + (_agent[1] ? 1 : 0);
}

// A micro-compartment is occupied if it is caseated or has at least one
// agent or has an external bacteria count > 0.
inline bool GridCell::isOccupied() const
{
	bool caseated = isCaseated();
	int agentCount = getNumberOfAgents();
	float extMtb = getExtMtb();

	bool result = caseated || (agentCount > 0) || (extMtb > 0.0);
	return result;
}


#endif /* GRIDCELL_H */
