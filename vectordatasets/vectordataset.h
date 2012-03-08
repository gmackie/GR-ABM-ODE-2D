/*
 * vectordataset.h
 *
 *  Created on: 12-sep-2008
 *      Author: s030858
 */

#ifndef VECTORDATASET_H_
#define VECTORDATASET_H_

#include "simulation.h"
#include "vector.h"
#include "datasets/dataset.h"

class VectorDataset : public Dataset
{

public:
	VectorDataset();
	virtual ~VectorDataset();
	virtual void getVector(const Simulation* pSimulation, int row, int col, vec2f& res) const = 0;
	vec2f getVectorNN(const Simulation* pSimulation, float row, float col) const;
	vec2f getVectorNN(const Simulation* pSimulation, const vec2f& pos) const;
	virtual void getVectorNN(const Simulation* pSimulation, const vec2f& pos, vec2f& res) const;
	vec2f getVectorBL(const Simulation* pSimulation, float row, float col) const;
	vec2f getVectorBL(const Simulation* pSimulation, const vec2f& pos) const;
	virtual void getVectorBL(const Simulation* pSimulation, const vec2f& pos, vec2f& res) const;
};

inline VectorDataset::VectorDataset()
	: Dataset()
{
}

inline VectorDataset::~VectorDataset()
{
}

inline vec2f VectorDataset::getVectorNN(const Simulation* pSimulation, float x, float y) const
{
	vec2f pos;
	pos[0] = x;
	pos[1] = y;

	return getVectorNN(pSimulation, pos);
}

inline vec2f VectorDataset::getVectorNN(const Simulation* pSimulation, const vec2f& pos) const
{
	vec2f res;
	getVectorNN(pSimulation, pos, res);
	return res;
}

inline void VectorDataset::getVectorNN(const Simulation* pSimulation, const vec2f& pos, vec2f& res) const
{
	int col = (int)pos[0];
	int row = (int)pos[1];
  const Pos& dim = pSimulation->getSize();

	if (pos[0] - col >= 0.5)
	{
		col = moduloDIM(col + 1, dim.x);
	}
	if (pos[1] - row >= 0.5)
	{
		row = moduloDIM(row + 1, dim.y);
	}

	getVector(pSimulation, row, col, res);
}

inline vec2f VectorDataset::getVectorBL(const Simulation* pSimulation, float x, float y) const
{
	vec2f pos;
	pos[0] = x;
	pos[1] = y;

	return getVectorBL(pSimulation, pos);
}

inline vec2f VectorDataset::getVectorBL(const Simulation* pSimulation, const vec2f& pos) const
{
	vec2f res;
	getVectorBL(pSimulation, pos, res);
	return res;
}

inline void VectorDataset::getVectorBL(const Simulation* pSimulation, const vec2f& pos, vec2f& res) const
{
	int col = (int)pos[0];
	int row = (int)pos[1];

  const Pos& dim = pSimulation->getSize();

	// make sure that we do not index out of bounds
	int sucRow = moduloDIM(row + 1, dim.x);
	int sucCol = moduloDIM(col + 1, dim.y);

	// first get the four vectors of the cell
	vec2f p1;
	getVector(pSimulation, row, col, p1);
	vec2f p2;
	getVector(pSimulation, row, sucCol, p2);
	vec2f p3;
	getVector(pSimulation, sucRow, col, p3);
	vec2f p4;
	getVector(pSimulation, sucRow, sucCol, p4);

	// interpolate in the x-direction
	//float t = pos[0] - (int)pos[0];
	float t = pos[0] - col;

	vec2f q1 = p1 + t * (p2 - p1);
	vec2f q2 = p3 + t * (p4 - p3);

	// interpolate in the y-direction
	//t = pos[1] - (int)pos[1];
	t = pos[1] - row;

	res = q1 + t * (q2 - q1);
}

#endif /* VECTORDATASET_H_ */
