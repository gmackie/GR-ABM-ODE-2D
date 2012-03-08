/*
 * vectorgrid.h
 *
 *  Created on: 21-sep-2008
 *      Author: S030858
 */

#ifndef VECTORGRID_H_
#define VECTORGRID_H_

#include "grviz.h"
#include "vector.h"
#include "vectordataset.h"
#include "datasets/grid.h"
#include <fstream>

struct VectorGridItem
{
	vec2f origin;
	vec2f vector;
};

class VectorGrid : public Grid
{
private:
	std::vector<VectorGridItem> _grid;

public:
	VectorGrid(size_t dim);
	~VectorGrid();
	void evaluate(const Simulation* pSimulation, VectorDataset* pVectorDataset, bool useNN);
	int count() const;
	const std::vector<VectorGridItem>& getGrid() const;
	void serialize(std::ofstream& outFile) const;
};

inline void VectorGrid::evaluate(const Simulation* pSimulation, VectorDataset* pVectorDataset, bool useNN)
{
	assert(pVectorDataset);

	for (size_t i = 0; i < _grid.size(); i++)
	{
		VectorGridItem& item = _grid[i];
		if (useNN)
		{
			pVectorDataset->getVectorNN(pSimulation, item.origin, item.vector);
		}
		else
		{
			pVectorDataset->getVectorBL(pSimulation, item.origin, item.vector);
		}
	}
}

inline int VectorGrid::count() const
{
	return (int)_grid.size();
}

inline const std::vector<VectorGridItem>& VectorGrid::getGrid() const
{
	return _grid;
}

#endif /* VECTORGRID_H_ */
