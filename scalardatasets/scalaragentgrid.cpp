/*
 * scalaragentgrid.cpp
 *
 *  Created on: 25-nov-2009
 *      Author: S030858
 */

#include "scalaragentgrid.h"
#include "simulation.h"
#include <fstream>

ScalarAgentGrid::ScalarAgentGrid()
	: ScalarAgentGridBase()
	, _grid(_DIM * _DIM)
{
	for (int i = 0; i < _DIM; i++) // row
	{
		for (int j = 0; j < _DIM; j++) // col
		{
			int idx = j + i * _DIM;
			_grid[idx]._bitMask = 0;
			_grid[idx]._nKillings = 0;
			_grid[idx]._pAgent[0] = _grid[idx]._pAgent[1] = NULL;
			_grid[idx]._attractant = 0;
			_grid[idx]._TNF = 0;
			_grid[idx]._CCL2 = 0;
			_grid[idx]._CCL5 = 0;
			_grid[idx]._CXCL9 = 0;
		}
	}
}

ScalarAgentGrid::~ScalarAgentGrid()
{
}

void ScalarAgentGrid::evaluate(const Simulation* pSimulation)
{
	const GrGrid& grGrid = pSimulation->getGrGrid();
	_macList = pSimulation->getMacList();
	_tgamList = pSimulation->getTgamList();
	_tcytList = pSimulation->getTcytList();
	_tregList = pSimulation->getTregList();

	for (int i = 0; i < _DIM; i++)
	{
		for (int j = 0; j < _DIM; j++)
		{
			const GridCell& cell = grGrid(i, j);
			ScalarAgentItem& item = _grid[i * _DIM + j];

			item._bitMask = 0;
			item._nKillings = cell.getNrKillings();
			item._nRecruitments = cell.getNrRecruitments();
			item._nSecretions = cell.getNrSecretions();
			item._attractant = cell.getMacAttractant();
			item._TNF = cell.getTNF();
			item._CCL2 = cell.getCCL2();
			item._CCL5 = cell.getCCL5();
			item._CXCL9 = cell.getCXCL9();
			item._extMtb = cell.getExtMtb();
			item._pAgent[0] = NULL;
			item._pAgent[1] = NULL;

			if (cell.isCaseated())
			{
				item._bitMask |= SET_BIT(_bitCas);
			}

			if (cell.isSource())
			{
				item._bitMask |= SET_BIT(_bitSrc);

				if (GrSimulation::MacRecruitmentThreshold(&cell))
				{
					item._bitMask |= SET_BIT(_bitSrcMac);
				}
				if (GrSimulation::TgamRecruitmentThreshold(&cell))
				{
					item._bitMask |= SET_BIT(_bitSrcTgam);
				}
				if (GrSimulation::TcytRecruitmentThreshold(&cell))
				{
					item._bitMask |= SET_BIT(_bitSrcTcyt);
				}
				if (GrSimulation::TregRecruitmentThreshold(&cell))
				{
					item._bitMask |= SET_BIT(_bitSrcTreg);
				}
			}

			if (cell.getExtMtb() >= 1)
				item._bitMask |= SET_BIT(_bitExtMtb);
		}
	}

	for (MacList::const_iterator it = _macList.begin(); it != _macList.end(); it++)
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
		case MAC_RESTING:
			item._bitMask |= (SET_BIT(_bitMac) | SET_BIT(_bitMacResting));
			break;
		case MAC_INFECTED:
			item._bitMask |= (SET_BIT(_bitMac) | SET_BIT(_bitMacInfected));
			break;
		case MAC_CINFECTED:
			item._bitMask |= (SET_BIT(_bitMac) | SET_BIT(_bitMacCInfected));
			break;
		case MAC_ACTIVE:
			item._bitMask |= (SET_BIT(_bitMac) | SET_BIT(_bitMacActive));
			break;
		case MAC_DEAD:
		//	item._bitMask |= (SET_BIT(_bitMac) | SET_BIT(_bitDead));
			break;
		}
	}

	for (TgamList::const_iterator it = _tgamList.begin(); it != _tgamList.end(); it++)
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
		case TGAM_DOWN_REGULATED:
			item._bitMask |= (SET_BIT(_bitTgam) | SET_BIT(_bitTgamDownRegulated));
			break;
		case TGAM_ACTIVE:
			item._bitMask |= (SET_BIT(_bitTgam) | SET_BIT(_bitTgamActive));
			break;
		case TGAM_DEAD:
		//	item._bitMask |= (SET_BIT(_bitTgam) | SET_BIT(_bitDead));
			break;
		}
	}

	for (TcytList::const_iterator it = _tcytList.begin(); it != _tcytList.end(); it++)
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
		case TCYT_DOWN_REGULATED:
			item._bitMask |= (SET_BIT(_bitTcyt) | SET_BIT(_bitTcytDownRegulated));
			break;
		case TCYT_ACTIVE:
			item._bitMask |= (SET_BIT(_bitTcyt) | SET_BIT(_bitTcytActive));
			break;
		case TCYT_DEAD:
		//	item._bitMask| = (SET_BIT(_bitTcyt) | SET_BIT(_bitDead));
			break;
		}
	}

	for (TregList::const_iterator it = _tregList.begin(); it != _tregList.end(); it++)
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
		case TREG_ACTIVE:
			item._bitMask |= (SET_BIT(_bitTreg) | SET_BIT(_bitTregResting));
			break;
		case TREG_DEAD:
		//	item._bitMask |= (SET_BIT(_bitTreg) | SET_BIT(_bitDead));
			break;
		}
	}
}
