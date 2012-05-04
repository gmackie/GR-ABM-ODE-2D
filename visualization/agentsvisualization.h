/*
 * agentsvisualization.h
 *
 *  Created on: 25-nov-2009
 *      Author: El-Kebir, M.
 */

#ifndef AGENTSVISUALIZATION_H_
#define AGENTSVISUALIZATION_H_

#include "visualization.h"
#include "scalardatasets/scalaragentgrid.h"

class AgentsVisualization : public Visualization
{
private:
	const ScalarAgentGrid* _pScalarAgentGrid;
	bool _drawGrid;
	bool _drawMacResting;
	bool _drawMacInfected;
	bool _drawMacCInfected;
	bool _drawMacActive;
	bool _drawTgam;
	bool _drawTcyt;
	bool _drawTreg;
	bool _drawCas;
	bool _drawSrc;
	bool _drawExtMtb;
	float _gridAlpha;
	float _gridHeight;
	bool _drawSrcMac;
	bool _drawSrcTgam;
	bool _drawSrcTcyt;
	bool _drawSrcTreg;
        int _drawm1m2;
        double _m1m2thres;
	int _selRow;
	int _selCol;
	void drawGrid() const;
	void drawQuad(int i, int j) const;
	void drawCross(int i, int j) const;
	void drawMark(int i, int j) const;

public:
	AgentsVisualization(int DIM, const ScalarAgentGrid* pScalarAgentGrid);
	virtual ~AgentsVisualization();
	void visualize(bool blend, const Simulation* pSimulation, const ColorMap* pColorMap) const;
	bool getDrawGrid() const;
	void setDrawGrid(bool drawGrid);
	float getGridAlpha() const;
	void setGridAlpha(float alpha);
	void setDrawMacResting(bool value);
	void setDrawMacInfected(bool value);
	void setDrawMacCInfected(bool value);
	void setDrawMacActive(bool value);
	void setDrawTgam(bool value);
	void setDrawTcyt(bool value);
	void setDrawTreg(bool value);
	void setDrawCas(bool value);
	void setDrawSrc(bool value);
	void setDrawExtMtb(bool value);
	void setGridHeight(float gridHeight);
	void setPredicates(bool mac, bool tgam, bool tcyt, bool treg);
    void setSelection(int row, int col);
    void setDrawM1M2(int value, double threshold);
};

inline void AgentsVisualization::setSelection(int row, int col)
{
	if (_selRow == row && _selCol == col)
	{
		_selRow = _selCol = -1;
	}
	else
	{
		_selRow = row;
		_selCol = col;
	}
}

inline void AgentsVisualization::setPredicates(bool mac, bool tgam, bool tcyt, bool treg)
{
	_drawSrcMac = mac;
	_drawSrcTgam = tgam;
	_drawSrcTcyt = tcyt;
	_drawSrcTreg = treg;
}

inline void AgentsVisualization::setGridHeight(float gridHeight)
{
	_gridHeight = gridHeight;
}

inline void AgentsVisualization::setDrawMacResting(bool value)
{
	_drawMacResting = value;
}

inline void AgentsVisualization::setDrawMacInfected(bool value)
{
	_drawMacInfected = value;
}

inline void AgentsVisualization::setDrawMacCInfected(bool value)
{
	_drawMacCInfected = value;
}

inline void AgentsVisualization::setDrawMacActive(bool value)
{
	_drawMacActive = value;
}

inline void AgentsVisualization::setDrawTgam(bool value)
{
	_drawTgam = value;
}

inline void AgentsVisualization::setDrawTcyt(bool value)
{
	_drawTcyt = value;
}

inline void AgentsVisualization::setDrawTreg(bool value)
{
	_drawTreg = value;
}

inline void AgentsVisualization::setDrawCas(bool value)
{
	_drawCas = value;
}

inline void AgentsVisualization::setDrawSrc(bool value)
{
	_drawSrc = value;
}

inline void AgentsVisualization::setDrawExtMtb(bool value)
{
	_drawExtMtb = value;
}

inline void AgentsVisualization::drawCross(int i, int j) const
{
	/*glBegin(GL_LINES);
		glVertex3f((j + .5) * deltaX, i * deltaY, _gridHeight);
		glVertex3f((j + .5) * deltaX, (i + 1) * deltaY, _gridHeight);
		glVertex3f((j + .1) * deltaX, (i + .65) * deltaY, _gridHeight);
		glVertex3f((j + .9) * deltaX, (i + .65) * deltaY, _gridHeight);
	glEnd();*/
	glBegin(GL_LINES);
		glVertex3f(j * _deltaX, i * _deltaY, 0.0f);
		glVertex3f((j + 1) * _deltaX, (i + 1) * _deltaY, 0.0f);
		glVertex3f(j * _deltaX, (i + 1) * _deltaY, 0.0f);
		glVertex3f((j + 1) * _deltaX, i * _deltaY, 0.0f);
	glEnd();

}

inline void AgentsVisualization::drawMark(int i, int j) const
{
	glBegin(GL_LINES);
		glVertex3f((j + 0.5) * _deltaX, i * _deltaY, 0.0f);
		glVertex3f((j + 0.5) * _deltaX, (i + 1) * _deltaY, 0.0f);
		glVertex3f(j * _deltaX, (i + 0.5) * _deltaY, 0.0f);
		glVertex3f((j + 1) * _deltaX, (i + 0.5) * _deltaY, 0.0f);
	glEnd();
}

inline void AgentsVisualization::drawQuad(int i, int j) const
{
	glBegin(GL_QUADS);
		glVertex3f(j * _deltaX, i * _deltaY, _gridHeight);
		glVertex3f(j * _deltaX, (i + 1) * _deltaY, _gridHeight);
		glVertex3f((j + 1) * _deltaX, (i + 1) * _deltaY, _gridHeight);
		glVertex3f((j + 1) * _deltaX, i * _deltaY, _gridHeight);
	glEnd();
}

inline bool AgentsVisualization::getDrawGrid() const
{
	return _drawGrid;
}

inline void AgentsVisualization::setDrawGrid(bool drawGrid)
{
	_drawGrid = drawGrid;
}

inline float AgentsVisualization::getGridAlpha() const
{
	return _gridAlpha;
}

inline void AgentsVisualization::setGridAlpha(float alpha)
{
	_gridAlpha = alpha;
}
inline void AgentsVisualization::setDrawM1M2(int value, double threshold)
{
    _drawm1m2 = value;
    _m1m2thres = threshold;
}


#endif /* AGENTSVISUALIZATION_H_ */
