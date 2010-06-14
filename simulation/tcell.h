/*
 * tcell.h
 *
 *  Created on: 13-nov-2009
 *      Author: M. El-Kebir
 */

#ifndef TCELL_H
#define TCELL_H

#include "agent.h"

class Tcell : public Agent
{
public:
	Tcell(int birthTime, int row, int col);
	virtual ~Tcell();
	void moveTcell(GrGrid& grid, bool ccl2, bool ccl5, bool cxcl9);
	static bool isTcell(const Agent* pAgent);
	virtual void serialize(std::ostream& out) const;
	virtual void deserialize(std::istream& in);
};

inline bool Tcell::isTcell(const Agent* pAgent)
{
	return pAgent && pAgent->getAgentType() != MAC;
}

#endif /* TCELL_H */
