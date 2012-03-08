/*
 * scalaragentgridbase.h
 *
 *  Created on: Mar 29, 2010
 *      Author: mohammed
 */

#ifndef SCALARAGENTGRIDBASE_H_
#define SCALARAGENTGRIDBASE_H_

#include "simulation.h"
#include "datasets/grid.h"

class ScalarAgentGridBase : public Grid
{
public:
	ScalarAgentGridBase(size_t dim);
	virtual ~ScalarAgentGridBase();
	virtual void evaluate(const Simulation* pSimulation) = 0;
};

inline ScalarAgentGridBase::ScalarAgentGridBase(size_t dim)
	: Grid(dim)
{
}

inline ScalarAgentGridBase::~ScalarAgentGridBase()
{
}

#endif /* SCALARAGENTGRIDBASE_H_ */
