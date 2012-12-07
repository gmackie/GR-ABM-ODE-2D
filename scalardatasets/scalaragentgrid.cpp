/*
 * scalaragentgrid.cpp
 *
 *  Created on: 25-nov-2009
 *      Author: S030858
 */

#include "scalaragentgrid.h"
#include "scalardataset.h"
#include "simulation.h"
#include "simulation/recruitmentbase.h"
#include "simulation/recruitmentprob.h"
#include <fstream>

ScalarAgentGrid::ScalarAgentGrid(size_t _DIM)
  : ScalarAgentGridBase(_DIM)
  , _grid(_DIM * _DIM)
{
}

ScalarAgentGrid::~ScalarAgentGrid()
{
}

template<typename T>
static void copy_list(const std::vector<T*>& src, std::vector<T>& dest)
{
  dest.reserve(src.size()); // Reduce the number of allocations we need to do
  dest.clear();
  for(size_t i=0;i<src.size();i++)
    dest.push_back(*(src[i]));
}

void ScalarAgentGrid::evaluate(const Simulation* pSimulation)
{
  copy_list(pSimulation->getMacList(), _macList);
  copy_list(pSimulation->getTgamList(), _tgamList);
  copy_list(pSimulation->getTcytList(), _tcytList);
  copy_list(pSimulation->getTregList(), _tregList);

  const GrGrid& grid = pSimulation->getGrGrid();

  for (int i = 0; i < _DIM; i++)
    {
      for (int j = 0; j < _DIM; j++)
        {
          const Pos p(i, j);
          ScalarAgentItem& item = _grid[i * _DIM + j];

          item._bitMask = 0;
          item._nKillings = grid.nKillings(p);
          item._nRecruitments = grid.nRecruitments(p);
          item._nRecruitmentsMac = grid.nRecruitmentsMac(p);
          item._nRecruitmentsTgam = grid.nRecruitmentsTgam(p);
          item._nRecruitmentsTcyt = grid.nRecruitmentsTcyt(p);
          item._nRecruitmentsTreg = grid.nRecruitmentsTreg(p);
          item._nSecretions = grid.nSecretions(p);
          item._attractant = grid.macAttractant(p);
          item._TNF = grid.TNF(p);
          item._CCL2 = grid.CCL2(p);
          item._CCL5 = grid.CCL5(p);
          item._CXCL9 = grid.CXCL9(p);
          item._shedTNFR2 = grid.shedTNFR2(p);
          item._il10 = grid.il10(p);
          item._extMtb = grid.extMTB(p);
          item._growthRate = grid.growthRate(p);
          item._pAgent[0] = NULL;
          item._pAgent[1] = NULL;

          if (grid.isCaseated(p))
            {
              item._bitMask |= SET_BIT(_bitCas);
            }

          if (grid.isSource(p))
            {
              item._bitMask |= SET_BIT(_bitSrc);

              if (RecruitmentProb::MacThresholdRecNew(grid, p))
                {
                  item._bitMask |= SET_BIT(_bitSrcMac);
                }
              if (RecruitmentProb::TgamThresholdRecNew(grid, p))
                {
                  item._bitMask |= SET_BIT(_bitSrcTgam);
                }
              if (RecruitmentProb::TcytThresholdRecNew(grid, p))
                {
                  item._bitMask |= SET_BIT(_bitSrcTcyt);
                }
              if (RecruitmentProb::TregThresholdRecNew(grid, p))
                {
                  item._bitMask |= SET_BIT(_bitSrcTreg);
                }
            }

          if (grid.extMTB(p) >= 1)
            item._bitMask |= SET_BIT(_bitExtMtb);
        }
    }

  for (typeof(_macList.begin()) it = _macList.begin(); it != _macList.end(); it++)
    {
      const Mac& mac = *it;
      ScalarAgentItem& item = _grid[mac.getRow() * _DIM + mac.getCol()];

      if (!item._pAgent[0])
        {
          item._pAgent[0] = &mac;
        }
      else
        {
          assert(!item._pAgent[1]);
          item._pAgent[1] = &mac;
        }

      if (mac.getNFkB())
        {
          item._bitMask |= SET_BIT(_bitMacNFkB);
        }
      if (mac.getStat1())
        {
          item._bitMask |= SET_BIT(_bitMacStat1);
        }

      item._bitMask |= SET_BIT(_bitMac);
      switch (mac.getState())
        {
        case Mac::MAC_RESTING:
          item._bitMask |= (SET_BIT(_bitMacResting));
          break;
        case Mac::MAC_INFECTED:
          item._bitMask |= (SET_BIT(_bitMacInfected));
          break;
        case Mac::MAC_CINFECTED:
          item._bitMask |= (SET_BIT(_bitMacCInfected));
          break;
        case Mac::MAC_ACTIVE:
          item._bitMask |= (SET_BIT(_bitMacActive));
          break;
        case Mac::MAC_DEAD:
          item._bitMask &= (~SET_BIT(_bitMac));
          break;
        }
    }

  for (typeof(_tgamList.begin()) it = _tgamList.begin(); it != _tgamList.end(); it++)
    {
      const Tgam& tgam = *it;
      ScalarAgentItem& item = _grid[tgam.getRow() * _DIM + tgam.getCol()];

      if (!item._pAgent[0])
        {
          item._pAgent[0] = &tgam;
        }
      else
        {
          assert(!item._pAgent[1]);
          item._pAgent[1] = &tgam;
        }

      item._bitMask |= SET_BIT(_bitTgam);
      switch (tgam.getState())
        {
        case Tgam::TGAM_DOWN_REGULATED:
          item._bitMask |= (SET_BIT(_bitTgamDownRegulated));
          break;
        case Tgam::TGAM_ACTIVE:
          item._bitMask |= (SET_BIT(_bitTgamActive));
          break;
        case Tgam::TGAM_DEAD:
          item._bitMask &= ~SET_BIT(_bitTgam);
          break;
        }
    }

  for (typeof(_tcytList.begin()) it = _tcytList.begin(); it != _tcytList.end(); it++)
    {
      const Tcyt& tcyt = *it;
      ScalarAgentItem& item = _grid[tcyt.getRow() * _DIM + tcyt.getCol()];

      if (!item._pAgent[0])
        {
          item._pAgent[0] = &tcyt;
        }
      else
        {
          assert(!item._pAgent[1]);
          item._pAgent[1] = &tcyt;
        }

      item._bitMask |= SET_BIT(_bitTcyt);
      switch (tcyt.getState())
        {
        case Tcyt::TCYT_DOWN_REGULATED:
          item._bitMask |= ( SET_BIT(_bitTcytDownRegulated));
          break;
        case Tcyt::TCYT_ACTIVE:
          item._bitMask |= ( SET_BIT(_bitTcytActive));
          break;
        case Tcyt::TCYT_DEAD:
            item._bitMask &= (~SET_BIT(_bitTcyt));
          break;
        }
    }

  for (typeof(_tregList.begin()) it = _tregList.begin(); it != _tregList.end(); it++)
    {
      const Treg& treg = *it;
      ScalarAgentItem& item = _grid[treg.getRow() * _DIM + treg.getCol()];

      if (!item._pAgent[0])
        {
          item._pAgent[0] = &treg;
        }
      else
        {
          assert(!item._pAgent[1]);
          item._pAgent[1] = &treg;
        }

      item._bitMask |= (SET_BIT(_bitTreg));
      switch (treg.getState())
        {
        case Treg::TREG_ACTIVE:
          item._bitMask |= (SET_BIT(_bitTreg) | SET_BIT(_bitTregResting));
          break;
        case Treg::TREG_DEAD:
          item._bitMask &= (~SET_BIT(_bitTreg));
          break;
        }
    }
}
