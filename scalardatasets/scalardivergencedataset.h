/*
 * scalardivergencedataset.h
 *
 *  Created on: 20-okt-2008
 *      Author: s030858
 */

#ifndef SCALARDIVERGENCEDATASET_H_
#define SCALARDIVERGENCEDATASET_H_

#include "scalardataset.h"
#include "vectordatasets/vectordataset.h"
#include "vectordatasets/vector.h"

class ScalarDivergenceDataset : public ScalarDataset
{
private:
	VectorDataset* _pVectorDataset;

public:
	ScalarDivergenceDataset(VectorDataset* pVectorDataset);
	virtual ~ScalarDivergenceDataset();
	virtual float getScalar(const Simulation* pSimulation, int row, int col) const;
	virtual float getScalarNN(const Simulation* pSimulation, float x, float y) const;
	virtual float getScalarBL(const Simulation* pSimulation, float x, float y) const;
};

inline ScalarDivergenceDataset::ScalarDivergenceDataset(VectorDataset* pVectorDataset)
	: ScalarDataset()
	, _pVectorDataset(pVectorDataset)
{
	assert(pVectorDataset);
}

inline ScalarDivergenceDataset::~ScalarDivergenceDataset()
{
	delete _pVectorDataset;
}

inline float ScalarDivergenceDataset::getScalar(const Simulation* pSimulation, int row, int col) const
{
	row = moduloDIM(pSimulation->getSize().y, row);
	col = moduloDIM(pSimulation->getSize().x, col);

	vec2f vectorLeft, vectorRight;
	_pVectorDataset->getVector(pSimulation, row, moduloDIM(pSimulation->getSize().x, col - 1), vectorLeft);
	_pVectorDataset->getVector(pSimulation, row, moduloDIM(pSimulation->getSize().x, col + 1), vectorRight);

	vec2f vectorDown, vectorUp;
	_pVectorDataset->getVector(pSimulation, moduloDIM(pSimulation->getSize().y, row - 1), col, vectorDown);
	_pVectorDataset->getVector(pSimulation, moduloDIM(pSimulation->getSize().y, row + 1), col, vectorUp);

	float dX = 0.5f * (vectorRight[0] - vectorLeft[0]);
	float dY = 0.5f * (vectorUp[1] - vectorDown[1]);

	return dX + dY;
}

inline float ScalarDivergenceDataset::getScalarNN(const Simulation* pSimulation, float x, float y) const
{
	float deltaX = 1;
	float deltaY = 1;

	vec2f vectorLeft, vectorRight;
	vectorLeft = _pVectorDataset->getVectorNN(pSimulation, moduloDIM(pSimulation->getSize().x, x - deltaX < _MIN_X ? _MAX_X : x - deltaX), y);
	vectorRight = _pVectorDataset->getVectorNN(pSimulation, moduloDIM(pSimulation->getSize().x, x + deltaX > _MAX_X ? _MIN_X : x + deltaX), y);

	vec2f vectorDown, vectorUp;
	vectorDown = _pVectorDataset->getVectorNN(pSimulation, x, moduloDIM(pSimulation->getSize().y, y - deltaY < _MIN_Y ? _MAX_Y : y - deltaY));
	vectorUp = _pVectorDataset->getVectorNN(pSimulation, x, moduloDIM(pSimulation->getSize().y, y + deltaY > _MAX_Y ? _MIN_Y : y + deltaY));

	float dX = (vectorRight[0] - vectorLeft[0]) / (2.0f * deltaX);
	float dY = (vectorUp[1] - vectorDown[1]) / (2.0f * deltaY);

	return dX + dY;
}

inline float ScalarDivergenceDataset::getScalarBL(const Simulation* pSimulation, float x, float y) const
{
	float deltaX = 1;
	float deltaY = 1;

	vec2f vectorLeft, vectorRight;
	vectorLeft = _pVectorDataset->getVectorBL(pSimulation, moduloDIM(pSimulation->getSize().x, x - deltaX < _MIN_X ? _MAX_X : x - deltaX), y);
	vectorRight = _pVectorDataset->getVectorBL(pSimulation, moduloDIM(pSimulation->getSize().x, x + deltaX > _MAX_X ? _MIN_X : x + deltaX), y);

	vec2f vectorDown, vectorUp;
	vectorDown = _pVectorDataset->getVectorBL(pSimulation, x, moduloDIM(pSimulation->getSize().y, y - deltaY < _MIN_Y ? _MAX_Y : y - deltaY));
	vectorUp = _pVectorDataset->getVectorBL(pSimulation, x, moduloDIM(pSimulation->getSize().y, y + deltaY > _MAX_Y ? _MIN_Y : y + deltaY));

	float dX = (vectorRight[0] - vectorLeft[0]) / (2.0f * deltaX);
	float dY = (vectorUp[1] - vectorDown[1]) / (2.0f * deltaY);

	/*std::cout << "(" << x << "," << y << "):\t"
		<< vectorRight[0] << "\t" << vectorLeft[0] << "\t"
		<< vectorUp[1] << "\t" << vectorDown[1] << std::endl;*/

	return dX + dY;
}

#endif /* SCALARDIVERGENCEDATASET_H_ */
