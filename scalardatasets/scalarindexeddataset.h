#ifndef SCALARINDEXEDDATASET_H
#define SCALARINDEXEDDATASET_H
#include "scalardataset.h"

class ScalarIndexedDataset : public ScalarDataset {
private:
  GrGrid::GRID_IDX idx;
public:
  ScalarIndexedDataset(GrGrid::GRID_IDX i) : ScalarDataset(), idx(i) {}
  float getScalar(const Simulation* pSimulation, int row, int col) const;
};

float ScalarIndexedDataset::getScalar(const Simulation* pSimulation, int row, int col) const {
  static float ret=0;
  pSimulation->getGrGrid().getIndexedValue(idx, row, col, ret);
  return ret;
}

#endif // SCALARINDEXEDDATASET_H
