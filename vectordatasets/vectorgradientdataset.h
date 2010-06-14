/*
 * VectorGradientDataset.h
 *
 *  Created on: 6-nov-2008
 *      Author: s030858
 */

#ifndef VECTORGRADIENTDATASET_H_
#define VECTORGRADIENTDATASET_H_

#include "vector.h"
#include "scalardatasets/scalardataset.h"
#include "vectordataset.h"

class VectorGradientDataset : public VectorDataset
{
private:
	ScalarDataset* _pScalarDataset;
	bool _shallow;

public:
	VectorGradientDataset(ScalarDataset* pScalarDataset, bool shallow = false);
	virtual void getVector(const Simulation* pSimulation, int row, int col, vec2f& res) const;
	virtual void getVectorNN(const Simulation* pSimulation, const vec2f& pos, vec2f& res) const;
	virtual void getVectorBL(const Simulation* pSimulation, const vec2f& pos, vec2f& res) const;
	virtual ~VectorGradientDataset();
};

inline VectorGradientDataset::VectorGradientDataset(ScalarDataset* pScalarDataset, bool shallow)
	: VectorDataset()
	, _pScalarDataset(pScalarDataset)
	, _shallow(shallow)
{
}

inline VectorGradientDataset::~VectorGradientDataset()
{
	if (!_shallow)
		delete _pScalarDataset;
}

inline void VectorGradientDataset::getVector(const Simulation* pSimulation, int row, int col, vec2f& res) const
{
	row = moduloDIM(row);
	col = moduloDIM(col);

	float scalarLeft = _pScalarDataset->getScalar(pSimulation, row, moduloDIM(col - 1));
	float scalarRight = _pScalarDataset->getScalar(pSimulation, row, moduloDIM(col + 1));

	float scalarDown = _pScalarDataset->getScalar(pSimulation, moduloDIM(row - 1), col);
	float scalarUp = _pScalarDataset->getScalar(pSimulation, moduloDIM(row + 1), col);

	res[0] = 0.5f * (scalarRight - scalarLeft);
	res[1] = 0.5f * (scalarUp - scalarDown);
}

inline void VectorGradientDataset::getVectorNN(const Simulation* pSimulation, const vec2f& pos, vec2f& res) const
{
	float scalarLeft = _pScalarDataset->getScalarNN(pSimulation, moduloDIM(pos[0] - 1), pos[1]);
	float scalarRight = _pScalarDataset->getScalarNN(pSimulation, moduloDIM(pos[0] + 1), pos[1]);

	float scalarDown = _pScalarDataset->getScalarNN(pSimulation, pos[0], moduloDIM(pos[1] - 1));
	float scalarUp = _pScalarDataset->getScalarNN(pSimulation, pos[0], moduloDIM(pos[1] + 1));

	res[0] = (scalarRight - scalarLeft) / 2.0f;
	res[1] = (scalarUp - scalarDown) / 2.0f;
}

inline void VectorGradientDataset::getVectorBL(const Simulation* pSimulation, const vec2f& pos, vec2f& res) const
{
	float scalarLeft = _pScalarDataset->getScalarBL(pSimulation, moduloDIM(pos[0] - 1), pos[1]);
	float scalarRight = _pScalarDataset->getScalarBL(pSimulation, moduloDIM(pos[0] + 1), pos[1]);

	float scalarDown = _pScalarDataset->getScalarBL(pSimulation, pos[0], moduloDIM(pos[1] -	1));
	float scalarUp = _pScalarDataset->getScalarBL(pSimulation, pos[0], moduloDIM(pos[1] + 1));

	res[0] = (scalarRight - scalarLeft) / 2.0f;
	res[1] = (scalarUp - scalarDown) / 2.0f;
}

#endif /* VECTORGRADIENTDATASET_H_ */
