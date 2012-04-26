/*
 * recruitmentprob.cpp
 *
 *  Created on: Jun 21, 2010
 *      Author: mohammed
 */

#include "recruitmentprob.h"
#include "stat.h"

RecruitmentProb::RecruitmentProb()
{
}

RecruitmentProb::~RecruitmentProb()
{
}

void RecruitmentProb::recruit(GrSimulation& sim)
{
	//const int timeTcellRecEnabled = _PARAM(PARAM_TCELL_TIME_RECRUITMENT_ENABLED);
  GrGrid& grid = sim.getGrid();
	const std::vector<Pos>& sources = grid.getSources();
	Stats& stats = sim.getStats();

	for (PosVector::const_iterator it = sources.begin();
		it != sources.end();
		it++)
	{
		const Pos& pSource = *it;

		// if the source is caseated continue
		if (grid.isCaseated(pSource))
			continue;

		// update stats
		if (MacRecruitmentThreshold(grid, pSource))
			++stats.getNrSources<Mac>();
		if (TgamRecruitmentThreshold(grid, pSource))
			++stats.getNrSources<Tgam>();
		if (TcytRecruitmentThreshold(grid, pSource))
			++stats.getNrSources<Tcyt>();
		if (TregRecruitmentThreshold(grid, pSource))
			++stats.getNrSources<Treg>();

		// Randomly choose the order of mac and T cell recruitment, so there isn't a bias in favor of one type of cell.
		if (g_Rand.getReal() < 0.5)
		{
			// macrophage recruitment
			if (grid.getNumberOfAgents(pSource) < 2 && !grid.hasAgentType(MAC, pSource))
				recruitMac(sim, pSource);

			// T cell recruitment
			if (sim.getTCellRecruitmentBegun() && grid.getNumberOfAgents(pSource) < 2)
				recruitTcell(sim, pSource);
		}
		else
		{
			// T cell recruitment
			if (sim.getTCellRecruitmentBegun() && grid.getNumberOfAgents(pSource) < 2)
				recruitTcell(sim, pSource);

			// macrophage recruitment
			if (grid.getNumberOfAgents(pSource) < 2 && !grid.hasAgentType(MAC, pSource))
				recruitMac(sim, pSource);
		}
	}
}

void RecruitmentProb::recruitMac(GrSimulation& sim, const Pos& pSource)
{
  GrGrid& grid = sim.getGrid();
	assert(!grid.isCaseated(pSource) && grid.getNumberOfAgents(pSource) < 2 && !grid.hasAgentType(MAC, pSource));

	// if the number of macrophages on the grid is less than _INITIAL_NUMBER_OF_MACROPHAGES,
	// recruit a resting macrophage
	if (sim.getStats().getNrOfMacs() < _PARAM(PARAM_MAC_INIT_NUMBER))
	{
		Mac* newMac = sim.createMac(pSource.x, pSource.y,
			sim.getTime() - g_Rand.getInt(_PARAM(PARAM_MAC_AGE)), Mac::MAC_RESTING, false, false);
		if (sim.getNfkbDynamics())
		{
			// initialize NF-kB signaling from steady-state
			for (int i = 0; i < 21600; ++i)
				newMac->solveNFkBODEsEquilibrium(2);
		}
	}
	else
	{
		bool macThreshold = MacRecruitmentThreshold(grid, pSource);

		if (macThreshold && g_Rand.getReal() < _PARAM(PARAM_MAC_PROB_RECRUITMENT))
		{
      ++grid.nRecruitments(pSource);
			Mac* newMac = sim.createMac(pSource.x, pSource.y,
				sim.getTime() - g_Rand.getInt(_PARAM(PARAM_MAC_AGE)), Mac::MAC_RESTING, false, false);
			if (sim.getNfkbDynamics())
			{
				// initialize NF-kB signaling from steady-state
				for (int i = 0; i < 21600; ++i)
					newMac->solveNFkBODEsEquilibrium(2);
			}
		}
	}
}

void RecruitmentProb::recruitTcell(GrSimulation& sim, const Pos& pSource)
{
  GrGrid& grid = sim.getGrid();
	assert(!grid.isCaseated(pSource) && grid.getNumberOfAgents(pSource) < 2);

	bool tgamThreshold = TgamRecruitmentThreshold(grid, pSource);
	bool tcytThreshold = TcytRecruitmentThreshold(grid, pSource);
	bool tregThreshold = TregRecruitmentThreshold(grid, pSource);

	if (g_Rand.getReal() < _PARAM(PARAM_TCELL_PROB_RECRUITMENT))
	{
		double r = g_Rand.getReal();
		if (r < _PARAM(PARAM_TGAM_PROB_RECRUITMENT))
		{
			// recruit a Tgam cell if allowed
			if (tgamThreshold)
			{
				++grid.nRecruitments(pSource);
				sim.createTgam(pSource.x, pSource.y,
					sim.getTime() - g_Rand.getInt(_PARAM(PARAM_TCELL_AGE), 1), Tgam::TGAM_ACTIVE);
			}
		}
		else if (r < (_PARAM(PARAM_TGAM_PROB_RECRUITMENT) + _PARAM(PARAM_TCYT_PROB_RECRUITMENT)))
		{
			// recruit a Tcyt cell if allowed
			if (tcytThreshold)
			{
				++grid.nRecruitments(pSource);
				sim.createTcyt(pSource.x, pSource.y,
					sim.getTime() - g_Rand.getInt(_PARAM(PARAM_TCELL_AGE), 1), Tcyt::TCYT_ACTIVE);
			}
		}
		else
		{
			// recruit a Treg cell if allowed
			if (tregThreshold)
			{
				++grid.nRecruitments(pSource);
				sim.createTreg(pSource.x, pSource.y,
					sim.getTime() - g_Rand.getInt(_PARAM(PARAM_TCELL_AGE), 1), Treg::TREG_ACTIVE);
			}
		}
	}
}
