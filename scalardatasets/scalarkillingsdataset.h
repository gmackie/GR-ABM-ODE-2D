#ifndef SCALARKILLINGSDATASET_H
#define SCALARKILLINGSDATASET_H

#include "scalardataset.h"

class ScalarKillingDataset:public ScalarDataset
{
public:
    ScalarKillingDataset();
    virtual ~ScalarKillingDataset();
    virtual float getScalar(const Simulation *pSimulation, int row, int col) const;
};

inline ScalarKillingDataset::ScalarKillingDataset()
{
}

inline ScalarKillingDataset::~ScalarKillingDataset()
{
}

inline float ScalarKillingDataset::getScalar(const Simulation *pSimulation, int row, int col) const
{
    return pSimulation->getGrGrid().nKillings(row, col);
}

#endif // SCALARKILLINGSDATASET_H
