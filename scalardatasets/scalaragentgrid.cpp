/*
 * scalaragentgrid.cpp
 *
 *  Created on: 25-nov-2009
 *      Author: S030858
 */

#include "scalaragentgrid.h"
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

void ScalarAgentGrid::evaluate(const Simulation* pSimulation)
{
	_macList = pSimulation->getMacList();
	_tgamList = pSimulation->getTgamList();
	_tcytList = pSimulation->getTcytList();
	_tregList = pSimulation->getTregList();

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
			item._nSecretions = grid.nSecretions(p);
			item._attractant = grid.macAttractant(p);
			item._TNF = grid.TNF(p);
			item._CCL2 = grid.CCL2(p);
			item._CCL5 = grid.CCL5(p);
			item._CXCL9 = grid.CXCL9(p);
			item._shedTNFR2 = grid.shedTNFR2(p);
			item._il10 = grid.il10(p);
			item._extMtb = grid.extMTB(p);
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
		ScalarAgentItem& item = _grid[it->getRow() * _DIM + it->getCol()];

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

		switch (mac.getState())
		{
		case Mac::MAC_RESTING:
			item._bitMask |= (SET_BIT(_bitMac) | SET_BIT(_bitMacResting));
			break;
		case Mac::MAC_INFECTED:
			item._bitMask |= (SET_BIT(_bitMac) | SET_BIT(_bitMacInfected));
			break;
		case Mac::MAC_CINFECTED:
			item._bitMask |= (SET_BIT(_bitMac) | SET_BIT(_bitMacCInfected));
			break;
		case Mac::MAC_ACTIVE:
			item._bitMask |= (SET_BIT(_bitMac) | SET_BIT(_bitMacActive));
			break;
		case Mac::MAC_DEAD:
		//	item._bitMask |= (SET_BIT(_bitMac) | SET_BIT(_bitDead));
			break;
		}
	}

	for (typeof(_tgamList.begin()) it = _tgamList.begin(); it != _tgamList.end(); it++)
	{
		const Tgam& tgam = *it;
		ScalarAgentItem& item = _grid[it->getRow() * _DIM + it->getCol()];

		if (!item._pAgent[0])
		{
			item._pAgent[0] = &tgam;
		}
		else
		{
			assert(!item._pAgent[1]);
			item._pAgent[1] = &tgam;
		}

		switch (tgam.getState())
		{
		case Tgam::TGAM_DOWN_REGULATED:
			item._bitMask |= (SET_BIT(_bitTgam) | SET_BIT(_bitTgamDownRegulated));
			break;
		case Tgam::TGAM_ACTIVE:
			item._bitMask |= (SET_BIT(_bitTgam) | SET_BIT(_bitTgamActive));
			break;
		case Tgam::TGAM_DEAD:
		//	item._bitMask |= (SET_BIT(_bitTgam) | SET_BIT(_bitDead));
			break;
		}
	}

	for (typeof(_tcytList.begin()) it = _tcytList.begin(); it != _tcytList.end(); it++)
	{
		const Tcyt& tcyt = *it;
		ScalarAgentItem& item = _grid[it->getRow() * _DIM + it->getCol()];

		if (!item._pAgent[0])
		{
			item._pAgent[0] = &tcyt;
		}
		else
		{
			assert(!item._pAgent[1]);
			item._pAgent[1] = &tcyt;
		}

		switch (tcyt.getState())
		{
		case Tcyt::TCYT_DOWN_REGULATED:
			item._bitMask |= (SET_BIT(_bitTcyt) | SET_BIT(_bitTcytDownRegulated));
			break;
		case Tcyt::TCYT_ACTIVE:
			item._bitMask |= (SET_BIT(_bitTcyt) | SET_BIT(_bitTcytActive));
			break;
		case Tcyt::TCYT_DEAD:
		//	item._bitMask| = (SET_BIT(_bitTcyt) | SET_BIT(_bitDead));
			break;
		}
	}

	for (typeof(_tregList.begin()) it = _tregList.begin(); it != _tregList.end(); it++)
	{
		const Treg& treg = *it;
		ScalarAgentItem& item = _grid[it->getRow() * _DIM + it->getCol()];

		if (!item._pAgent[0])
		{
			item._pAgent[0] = &treg;
		}
		else
		{
			assert(!item._pAgent[1]);
			item._pAgent[1] = &treg;
		}

		switch (treg.getState())
		{
		case Treg::TREG_ACTIVE:
			item._bitMask |= (SET_BIT(_bitTreg) | SET_BIT(_bitTregResting));
			break;
		case Treg::TREG_DEAD:
		//	item._bitMask |= (SET_BIT(_bitTreg) | SET_BIT(_bitDead));
			break;
		}
	}
}
