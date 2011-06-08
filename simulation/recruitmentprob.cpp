/*
 * recruitmentprob.cpp
 *
 *  Created on: Jun 21, 2010
 *      Author: mohammed
 */

#include "recruitmentprob.h"
#include "grstat.h"

RecruitmentProb::RecruitmentProb()
{
}

RecruitmentProb::~RecruitmentProb()
{
}

void RecruitmentProb::recruit(GrSimulation& sim)
{
	const int timeTcellRecEnabled = _PARAM(PARAM_TCELL_TIME_RECRUITMENT_ENABLED);
	GridCellPtrVector& sources = sim.getGrid().getSources();
	GrStat& stats = sim.getStats();

	for (GridCellPtrVector::iterator it = sources.begin();
		it != sources.end();
		it++)
	{
		GridCell* pSource = *it;

		// if the source is caseated continue
		if (pSource->isCaseated())
			continue;

		// update stats
		if (MacRecruitmentThreshold(pSource))
			stats.incNrSourcesMac();
		if (TgamRecruitmentThreshold(pSource))
			stats.incNrSourcesTgam();
		if (TcytRecruitmentThreshold(pSource))
			stats.incNrSourcesTcyt();
		if (TregRecruitmentThreshold(pSource))
			stats.incNrSourcesTreg();

		// Randomly choose the order of mac and T cell recruitment, so there isn't a bias in favor of one type of cell.
		if (g_Rand.getReal() < 0.5)
		{
			// macrophage recruitment
			if (pSource->getNumberOfAgents() < 2 && !pSource->hasMac())
				recruitMac(sim, pSource);

			// T cell recruitment
			if (pSource->getNumberOfAgents() < 2 && sim.getTime() >= timeTcellRecEnabled)
				recruitTcell(sim, pSource);
		}
		else
		{
			// T cell recruitment
			if (pSource->getNumberOfAgents() < 2 && sim.getTime() >= timeTcellRecEnabled)
				recruitTcell(sim, pSource);

			// macrophage recruitment
			if (pSource->getNumberOfAgents() < 2 && !pSource->hasMac())
				recruitMac(sim, pSource);
		}
	}
}

void RecruitmentProb::recruitMac(GrSimulation& sim, GridCell* pSource)
{
	assert(!pSource->isCaseated() && pSource->getNumberOfAgents() < 2 && !pSource->hasMac());

	// if the number of macrophages on the grid is less than _INITIAL_NUMBER_OF_MACROPHAGES,
	// recruit a resting macrophage
	if (sim.getStats().getNrOfMac() < _PARAM(PARAM_MAC_INIT_NUMBER))
	{
		Mac* newMac = sim.createMac(pSource->getRow(), pSource->getCol(),
			sim.getTime() - g_Rand.getInt(_PARAM(PARAM_MAC_AGE)), MAC_RESTING, false, false);
		if (sim.getNfkbDynamics())
		{
			// initialize NF-kB signaling from steady-state
			for (int i = 0; i < 21600; ++i)
				newMac->solveNFkBODEsEquilibrium(2);
		}
	}
	else
	{
		bool macThreshold = MacRecruitmentThreshold(pSource);

		if (macThreshold && g_Rand.getReal() < _PARAM(PARAM_MAC_PROB_RECRUITMENT))
		{
			pSource->incNrRecruitments();
			Mac* newMac = sim.createMac(pSource->getRow(), pSource->getCol(),
				sim.getTime() - g_Rand.getInt(_PARAM(PARAM_MAC_AGE)), MAC_RESTING, false, false);
			if (sim.getNfkbDynamics())
			{
				// initialize NF-kB signaling from steady-state
				for (int i = 0; i < 21600; ++i)
					newMac->solveNFkBODEsEquilibrium(2);
			}
		}
	}
}

void RecruitmentProb::recruitTcell(GrSimulation& sim, GridCell* pSource)
{
	assert(!pSource->isCaseated() && pSource->getNumberOfAgents() < 2);

	bool tgamThreshold = TgamRecruitmentThreshold(pSource);
	bool tcytThreshold = TcytRecruitmentThreshold(pSource);
	bool tregThreshold = TregRecruitmentThreshold(pSource);

	if (g_Rand.getReal() < _PARAM(PARAM_TCELL_PROB_RECRUITMENT))
	{
		double r = g_Rand.getReal();
		if (r < _PARAM(PARAM_TGAM_PROB_RECRUITMENT))
		{
			// recruit a Tgam cell if allowed
			if (tgamThreshold)
			{
				pSource->incNrRecruitments();
				sim.createTgam(pSource->getRow(), pSource->getCol(),
					sim.getTime() - g_Rand.getInt(_PARAM(PARAM_TCELL_AGE), 1), TGAM_ACTIVE);
			}
		}
		else if (r < (_PARAM(PARAM_TGAM_PROB_RECRUITMENT) + _PARAM(PARAM_TCYT_PROB_RECRUITMENT)))
		{
			// recruit a Tcyt cell if allowed
			if (tcytThreshold)
			{
				pSource->incNrRecruitments();
				sim.createTcyt(pSource->getRow(), pSource->getCol(),
					sim.getTime() - g_Rand.getInt(_PARAM(PARAM_TCELL_AGE), 1), TCYT_ACTIVE);
			}
		}
		else
		{
			// recruit a Treg cell if allowed
			if (tregThreshold)
			{
				pSource->incNrRecruitments();
				sim.createTreg(pSource->getRow(), pSource->getCol(),
					sim.getTime() - g_Rand.getInt(_PARAM(PARAM_TCELL_AGE), 1), TREG_ACTIVE);
			}
		}
	}
}
