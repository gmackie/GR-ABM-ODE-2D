/*
 * GrSimulationGrid.h
 *
 *  Created on: Nov 23, 2010
 *      Author: pwolberg
 */

#ifndef GRSIMULATIONGRID_H_
#define GRSIMULATIONGRID_H_

#include "grgrid.h"

// A class to encapsulate management of the simulation grid for
// grid swapping during diffusion.
// Other than grid swapping during diffusion, all operations
// an a GrSimulationGrid object refer to/operate on the current grid.
class GrSimulationGrid
{
public:
  GrSimulationGrid(const Pos& dim);
  GrSimulationGrid(const GrSimulationGrid& g);
  virtual ~GrSimulationGrid();

  const Pos& getRange() const
  {
    return getGrid().getRange();
  }
  size_t getSize() const
  {
    return getGrid().getSize();
  }
  const Pos getCenter() const
  {
    return  getGrid().getCenter();
  }

  const GrGrid& getCurrentGrid() const;
  GrGrid& getCurrentGrid();
  const GrGrid& getNextGrid() const;
  GrGrid& getNextGrid();

  void updateNextGrid();
  void swap();

  // For backward compatibility with existing code.
  // These all refer to/operate on the current grid.
  const GrGrid& getGrid() const;
  GrGrid& getGrid();
  const PosVector& getSources();
  void initSources();
  template<typename Archive>
  void serialize(Archive& ar, const unsigned int version);

private:
  GrGrid *_pCurrentGrid;
  GrGrid *_pNextGrid;
};


inline const GrGrid& GrSimulationGrid::getCurrentGrid() const
{
  return (*_pCurrentGrid);
}
inline GrGrid& GrSimulationGrid::getCurrentGrid()
{
  return (*_pCurrentGrid);
}
inline const GrGrid& GrSimulationGrid::getNextGrid() const
{
  return (*_pNextGrid);
}
inline GrGrid& GrSimulationGrid::getNextGrid()
{
  return (*_pNextGrid);
}

inline const GrGrid& GrSimulationGrid::getGrid() const
{
  return (*_pCurrentGrid);
}

inline GrGrid& GrSimulationGrid::getGrid()
{
  return (*_pCurrentGrid);
}

inline const PosVector& GrSimulationGrid::getSources()
{
  return _pCurrentGrid->getSources();
}

inline void GrSimulationGrid::initSources()
{
  _pCurrentGrid->initSources();
}

// Make sure that the next grid has the same state as the current grid.
// This includes agent pointers in addition to chemical concentrations
// in each grid compartment.
// This is necessary before using grid swapping at the start of the
// diffusion loop.
inline void GrSimulationGrid::updateNextGrid()
{
  *_pNextGrid = *_pCurrentGrid;
}

// Swap the current and next grids.
// Done after the next state has been stored in the next grid based
// on the current state in the current.
inline void GrSimulationGrid::swap()
{
  std::swap(_pCurrentGrid, _pNextGrid);
}

template<typename Archive>
void GrSimulationGrid::serialize(Archive& ar, const unsigned int /*version*/) {
  ar & boost::serialization::make_nvp("CurrentGrid", *_pCurrentGrid);
}

#endif /* GRSIMULATIONGRID_H_ */
