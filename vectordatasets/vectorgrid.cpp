/*
 * vectorgrid.cpp
 *
 *  Created on: 21-sep-2008
 *      Author: S030858
 */

#include "vectorgrid.h"
#include "vectordataset.h"
#include "simulation.h"

VectorGrid::VectorGrid(size_t _DIM)
  : Grid(_DIM)
  , _grid(_DIM*_DIM)
{
  for (size_t i = 0; i < _DIM; i++) // row
    {
      for (size_t j = 0; j < _DIM; j++) // col
        {
          VectorGridItem& item = _grid[j + i * _DIM];

          item.origin[0] = j;
          item.origin[1] = i;
        }
    }
}

void VectorGrid::serialize(std::ofstream& outFile) const
{
  outFile << _DIM << ";" << _DIM << std::endl;

  for (int i = 0; i < _DIM; i++) // row
    {
      for (int j = 0; j < _DIM; j++) // col
        {
          if (j != 0)
            {
              outFile << ";";
            }
          const VectorGridItem& item = _grid[i * _DIM + j];
          outFile << item.vector[0];
        }
      outFile << std::endl;
    }

  outFile << std::endl;

  for (int i = 0; i < _DIM; i++) // row
    {
      for (int j = 0; j < _DIM; j++) // col
        {
          if (j != 0)
            {
              outFile << ";";
            }
          const VectorGridItem& item = _grid[i * _DIM + j];
          outFile << item.vector[1];
        }
      outFile << std::endl;
    }

  outFile << std::endl;
}

VectorGrid::~VectorGrid()
{
}

void VectorGrid::evaluate(const Simulation* pSimulation, VectorDataset* pVectorDataset, bool useNN)
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
