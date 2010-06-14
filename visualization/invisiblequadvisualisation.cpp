/*
 * invisiblequadvisualisation.cpp
 *
 *  Created on: 7-dec-2008
 *      Author: s030858
 */

#include "invisiblequadvisualisation.h"

InvisibleQuadVisualisation::InvisibleQuadVisualisation(int DIM)
	: Visualization(DIM)
{
}

InvisibleQuadVisualisation::~InvisibleQuadVisualisation()
{
}

void InvisibleQuadVisualisation::visualize(bool, const Simulation*, const ColorMap*) const
{
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);
	glBegin(GL_QUADS);
		glVertex3f(_MIN_X, _MIN_Y, 0.0f);
		glVertex3f(_MAX_X, _MIN_Y, 0.0f);
		glVertex3f(_MAX_X, _MAX_Y, 0.0f);
		glVertex3f(_MIN_X, _MAX_Y, 0.0f);
	glEnd();
	glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
}
