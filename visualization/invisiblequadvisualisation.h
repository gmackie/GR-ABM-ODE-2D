/*
 * invisiblequadvisualisation.h
 *
 *  Created on: 7-dec-2008
 *      Author: s030858
 */

#ifndef INVISIBLEQUADVISUALISATION_H_
#define INVISIBLEQUADVISUALISATION_H_

#include "visualization.h"

class InvisibleQuadVisualisation: public Visualization
{
public:
	InvisibleQuadVisualisation(int DIM);
	virtual ~InvisibleQuadVisualisation();
    void visualize(bool blend, const Simulation* pSimulation, const ColorMap* pColorMap) const;
};

#endif /* INVISIBLEQUADVISUALISATION_H_ */
