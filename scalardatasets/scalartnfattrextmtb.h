/*
 * scalarTnfAttrExtMtb.h
 *
 * A scalar dataset for determining a granuloma boundary.
 * It uses the sum of TNF, the Macrophage attractant secreted by a caseated region
 * and the extra-cellular bacteria count.
 *
 * TNF: the usual determinant of a granuloma boundary.
 *
 * Attractant: To make sure any caseated region on the edge of a granuloma is inside
 *             the boundary. Because there is low or no TNF near a caseated region,
 *             a caseated region won't be included within a granuloma based on TNF alone.
 *
 *  Extra-cellular bacteria count: A region of dense bacteria in a granuloma, usually in the
 *                                 center, will have no TNF or attractant, so this region will
 *                                 be excluded if only using TNF and attractant. This creates
 *                                 a doughnut shaped granuloma with a hole in the center with
 *                                 a boundary line around the hole. Including bacteria count
 *                                 should solve this problem.
 *
 *  Created on: Apr 27, 2010
 *      Author: Paul Wolberg
 */

#ifndef SCALARTNFATTREXTMTB_H_
#define SCALARTNFATTREXTMTB_H_

#include "scalardataset.h"

class ScalarTnfAttrExtMtb: public ScalarDataset {
public:
	ScalarTnfAttrExtMtb();
	virtual ~ScalarTnfAttrExtMtb();
	virtual float getScalar(const Simulation* pSimulation, int row, int col) const;
};

inline ScalarTnfAttrExtMtb::ScalarTnfAttrExtMtb()
	: ScalarDataset()
{
}

inline ScalarTnfAttrExtMtb::~ScalarTnfAttrExtMtb()
{
}

inline float ScalarTnfAttrExtMtb::getScalar(const Simulation* pSimulation, int row, int col) const
{
	return pSimulation->getGrGrid().TNF(row, col)
           + pSimulation->getGrGrid().macAttractant(row, col)
           + pSimulation->getGrGrid().extMTB(row, col);
}

#endif /* SCALARTNFATTREXTMTB_H_ */













