/*
 * ScalarCellDensityDataset.h
 *
 * A scalar dataset for determining a granuloma boundary by cell density.
 *
 * This uses an idea by Mohammad Fallahi
 *
 *  Created on: Apr 29, 2010
 *      Author: Paul Wolberg
 */

#ifndef SCALARCELLDENSITYDATASET_H_
#define SCALARCELLDENSITYDATASET_H_

#include "scalardataset.h"

class ScalarCellDensityDataset: public ScalarDataset
{
public:
  ScalarCellDensityDataset();
  virtual ~ScalarCellDensityDataset();
  virtual float getScalar(const Simulation* pSimulation, int row, int col) const;
};

inline ScalarCellDensityDataset::ScalarCellDensityDataset()
  : ScalarDataset()
{
}

inline ScalarCellDensityDataset::~ScalarCellDensityDataset()
{
}

inline float ScalarCellDensityDataset::getScalar(const Simulation* pSimulation, int row, int col) const
{
  float cellDensity = 0.0; // The cell density in the Moore neighborhood of the micro-compartment at (row, col).
  float mcCount = 0.0;     // The number of cells in the Moore neighborhood that are occupied - are caseated or have cells.

  // Check the micro-compartment's Moore neighborhood (including the micro-compartment itself).
  for (int i = -1; i <= 1; i++)
    {
      for (int j = -1; j <= 1; j++)
        {
          if (pSimulation->getGrGrid().isOccupied(moduloDIM(pSimulation->getSize().y, row + i), moduloDIM(pSimulation->getSize().x, col + j)))
            {
              mcCount++;
            }
        }
    }

  cellDensity = mcCount/MOORE_COUNT;

  return cellDensity;
}


#endif /* SCALARCELLDENSITYDATASET_H_ */
