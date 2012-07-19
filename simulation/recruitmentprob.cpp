/*
 * recruitmentprob.cpp
 *
 *  Created on: Jun 21, 2010
 *      Author: mohammed
 */

#include "recruitmentprob.h"
#include "grsimulation.h"
#include "stat.h"

RecruitmentProb::RecruitmentProb()
{
}

RecruitmentProb::~RecruitmentProb()
{
}

void RecruitmentProb::recruit(GrSimulation &sim)
{
    GrGrid& grid = sim.getGrid();
    const std::vector<Pos>& sources = grid.getSources();
    //Stats& stats = sim.getStats();

//    std::cout << "Recruit Loop" << std::endl;

    for (PosVector::const_iterator it = sources.begin(); it != sources.end(); it++)
    {

        if (grid.getNumberOfAgents(*it) == (int)GrGrid::MAX_AGENTS_PER_CELL || grid.isCaseated(*it))
        {
//            std::cout << "No Recruitment" << std::endl;
            continue;
        }

        for(size_t i=0; i<(GrGrid::MAX_AGENTS_PER_CELL - grid.getNumberOfAgents(*it)); i++)
        {
//            std::cout << "RecruitStep: " << i << std::endl;
            recruit(sim, *it);
        }
    }
}



void RecruitmentProb::recruit(GrSimulation &sim, const Pos& pSource)
{
    GrGrid& grid = sim.getGrid();
    Stats& stats = sim.getStats();

//    std::cout << "Recruit Loop" << std::endl;

    double runningProbability = 0.0;
    const double randomNumber = g_Rand.getReal();
    const double evenProb = 0.25;
    const bool pmac = PossibleRecruitMac(sim, pSource)  && !grid.hasAgentType(MAC, pSource);
    const bool ptgam = PossibleRecruitTgam(sim, pSource) && sim.getTCellRecruitmentBegun();
    const bool ptcyt = PossibleRecruitTcyt(sim, pSource) && sim.getTCellRecruitmentBegun();
    const bool ptreg = PossibleRecruitTreg(sim, pSource) && sim.getTCellRecruitmentBegun();
//    std::cout << "Possible Mac: " << pmac  << "  Possible Tgam: " << ptgam << "  Possible Tcyt: " << ptcyt << "  Possible Treg: " << ptreg << std::endl;

    // Update Recruitment Stats
    stats.getNrSources<Mac>() += pmac;
    stats.getNrSources<Tgam>() += ptgam;
    stats.getNrSources<Tcyt>() += ptcyt;
    stats.getNrSources<Treg>() += ptreg;

    const size_t countPossible = pmac + ptgam + ptcyt + ptreg;

    if (pmac)
    {
        runningProbability += (evenProb * pmac) + (evenProb / countPossible) * (!ptgam +
                            !ptcyt + !ptreg);

        if ((randomNumber < runningProbability))
        {
            recruitCellMac(sim, pSource);
//            std::cout << "Recruit Mac and Break" << std::endl;
            return;
        }

    }

    if (ptgam)
    {
        runningProbability += (evenProb * ptgam) + (evenProb / countPossible) * (!pmac +
                            !ptcyt + !ptreg);

        if ((randomNumber < runningProbability))
        {
            recruitCellTgam(sim, pSource);
//            std::cout << "Recruit Tgam and Break" << std::endl;
            return;
        }

    }

    if (ptcyt)
    {
        runningProbability += (evenProb * ptcyt) + (evenProb / countPossible) * (!ptgam +
                            !pmac + !ptreg);

        if ((randomNumber < runningProbability))
        {
            recruitCellTcyt(sim, pSource);
//            std::cout << "Recruit Tcyt and Break" << std::endl;
            return;
        }

    }

    if (ptreg)
    {
        runningProbability += (evenProb * ptreg) + (evenProb / countPossible) * (!ptgam +
                            !ptcyt + !ptreg);

        if ((randomNumber < runningProbability))
        {
            recruitCellTreg(sim, pSource);
//            std::cout << "Recruit Treg and Break" << std::endl;
            return;
        }

    }

}


//void RecruitmentProb::recruit(GrSimulation& sim)
//{
//    //const int timeTcellRecEnabled = _PARAM(PARAM_TCELL_TIME_RECRUITMENT_ENABLED);
//  GrGrid& grid = sim.getGrid();
//    const std::vector<Pos>& sources = grid.getSources();
//    Stats& stats = sim.getStats();

//    for (PosVector::const_iterator it = sources.begin();
//        it != sources.end();
//        it++)
//    {
//        const Pos& pSource = *it;

//        // if the source is caseated continue
//        if (grid.isCaseated(pSource))
//            continue;

//        // update stats
//        if (MacRecruitmentThreshold(grid, pSource))
//            ++stats.getNrSources<Mac>();
//        if (TgamRecruitmentThreshold(grid, pSource))
//            ++stats.getNrSources<Tgam>();
//        if (TcytRecruitmentThreshold(grid, pSource))
//            ++stats.getNrSources<Tcyt>();
//        if (TregRecruitmentThreshold(grid, pSource))
//            ++stats.getNrSources<Treg>();

//        // Randomly choose the order of mac and T cell recruitment, so there isn't a bias in favor of one type of cell.
//        if (g_Rand.getReal() < 0.5)
//        {
//            // macrophage recruitment
//            if (grid.getNumberOfAgents(pSource) < 2 && !grid.hasAgentType(MAC, pSource))
//                recruitMac(sim, pSource);

//            // T cell recruitment
//            if (sim.getTCellRecruitmentBegun() && grid.getNumberOfAgents(pSource) < 2)
//                recruitTcell(sim, pSource);
//        }
//        else
//        {
//            // T cell recruitment
//            if (sim.getTCellRecruitmentBegun() && grid.getNumberOfAgents(pSource) < 2)
//                recruitTcell(sim, pSource);

//            // macrophage recruitment
//            if (grid.getNumberOfAgents(pSource) < 2 && !grid.hasAgentType(MAC, pSource))
//                recruitMac(sim, pSource);
//        }
//    }
//}

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

//        double macThresholdNewValue = 0.0;

//        bool macBoolNew = MacThresholdRecNew(grid, pSource, macThresholdNewValue);

//        std::cout << "Function Value: " << macThresholdNewValue << "  Bool Value: " << macBoolNew << std::endl;


//        if (macBoolNew && g_Rand.getReal() < (_PARAM(PARAM_MAC_MAX_RECRUITMENT) * macThresholdNewValue))
//        {
//            ++grid.nRecruitments(pSource);
//            Mac* newMac = sim.createMac(pSource.x, pSource.y,
//                sim.getTime() - g_Rand.getInt(_PARAM(PARAM_MAC_AGE)), Mac::MAC_RESTING, false, false);
//            if (sim.getNfkbDynamics())
//            {
//                // initialize NF-kB signaling from steady-state
//                for (int i = 0; i < 21600; ++i)
//                    newMac->solveNFkBODEsEquilibrium(2);
//            }

//        }

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


inline void RecruitmentProb::recruitCellMac(GrSimulation &sim, const Pos &pSource)
{
    GrGrid& grid = sim.getGrid();

    assert(!grid.isCaseated(pSource) && grid.getNumberOfAgents(pSource) < 2 && !grid.hasAgentType(MAC, pSource));

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

inline void RecruitmentProb::recruitCellTgam(GrSimulation &sim, const Pos &pSource)
{
    GrGrid& grid = sim.getGrid();

    assert(!grid.isCaseated(pSource) && grid.getNumberOfAgents(pSource) < 2);

    ++grid.nRecruitments(pSource);
    sim.createTgam(pSource.x, pSource.y,
        sim.getTime() - g_Rand.getInt(_PARAM(PARAM_TCELL_AGE), 1), Tgam::TGAM_ACTIVE);
}


inline void RecruitmentProb::recruitCellTcyt(GrSimulation &sim, const Pos &pSource)
{
    GrGrid& grid = sim.getGrid();

    assert(!grid.isCaseated(pSource) && grid.getNumberOfAgents(pSource) < 2);

    ++grid.nRecruitments(pSource);
    sim.createTcyt(pSource.x, pSource.y,
        sim.getTime() - g_Rand.getInt(_PARAM(PARAM_TCELL_AGE), 1), Tcyt::TCYT_ACTIVE);
}

inline void RecruitmentProb::recruitCellTreg(GrSimulation &sim, const Pos &pSource)
{
    GrGrid& grid = sim.getGrid();

    assert(!grid.isCaseated(pSource) && grid.getNumberOfAgents(pSource) < 2);

    ++grid.nRecruitments(pSource);
    sim.createTreg(pSource.x, pSource.y,
        sim.getTime() - g_Rand.getInt(_PARAM(PARAM_TCELL_AGE), 1), Treg::TREG_ACTIVE);
}


bool RecruitmentProb::PossibleRecruitMac(GrSimulation& sim, const Pos& pSource)
{
    GrGrid& grid = sim.getGrid();

    if (!grid.isCaseated(pSource) && grid.getNumberOfAgents(pSource) < 2 && !grid.hasAgentType(MAC, pSource))
    {
        // if the number of macrophages on the grid is less than _INITIAL_NUMBER_OF_MACROPHAGES,
        // recruit a resting macrophage
        if (sim.getStats().getNrOfMacs() < _PARAM(PARAM_MAC_INIT_NUMBER))
            return true;

        double macProbValue = 0.0;
        bool macBool = MacThresholdRecNew(grid, pSource, macProbValue);

//        std::cout << "Mac: " << pSource  << "  "  << macProbValue  << "  " << macBool << std::endl;

        return (macBool && g_Rand.getReal() < (_PARAM(PARAM_MAC_MAX_RECRUITMENT) * macProbValue));
    }

    else
    {
        return false;
    }
}

bool RecruitmentProb::PossibleRecruitTgam(GrSimulation& sim, const Pos& pSource)
{
    GrGrid& grid = sim.getGrid();

    if (!grid.isCaseated(pSource) && grid.getNumberOfAgents(pSource) < 2)
    {
        double TgamProbValue = 0.0;
        bool TgamBool = TgamThresholdRecNew(grid, pSource, TgamProbValue);

        if (grid.getNumberOfAgents(pSource) == 1)
        {
            if (grid.hasAgentType(MAC, pSource))
                TgamBool = g_Rand.getReal() < _PARAM(PARAM_TCELL_PROB_MOVE_TO_MAC);
            else
                TgamBool = g_Rand.getReal() < _PARAM(PARAM_TCELL_PROB_MOVE_TO_TCELL);
        }

//        if (pSource == Pos (53,48))
//            std::cout << "Tgam: " << pSource  << "  "  << TgamProbValue  << "  " << TgamBool << std::endl;


        return (TgamBool && g_Rand.getReal() < (_PARAM(PARAM_TGAM_MAX_RECRUITMENT) * TgamProbValue));

    }
    else
    {
        return false;
    }
}

bool RecruitmentProb::PossibleRecruitTcyt(GrSimulation& sim, const Pos& pSource)
{
    GrGrid& grid = sim.getGrid();

    if (!grid.isCaseated(pSource) && grid.getNumberOfAgents(pSource) < 2)
    {
        double TcytProbValue = 0.0;
        bool TcytBool = TcytThresholdRecNew(grid, pSource, TcytProbValue);


        if (grid.getNumberOfAgents(pSource) == 1)
        {
            if (grid.hasAgentType(MAC, pSource))
                TcytBool = g_Rand.getReal() < _PARAM(PARAM_TCELL_PROB_MOVE_TO_MAC);
            else
                TcytBool = g_Rand.getReal() < _PARAM(PARAM_TCELL_PROB_MOVE_TO_TCELL);
        }

//        if (pSource == Pos (53,48))
//            std::cout << "Tcyt: " << pSource  << "  "  << TcytProbValue << "  " << TcytBool << std::endl;

        return (TcytBool && g_Rand.getReal() < (_PARAM(PARAM_TCYT_MAX_RECRUITMENT) * TcytProbValue));
    }
    else
    {
        return false;
    }
}

bool RecruitmentProb::PossibleRecruitTreg(GrSimulation& sim, const Pos& pSource)
{
    GrGrid& grid = sim.getGrid();

    if (!grid.isCaseated(pSource) && grid.getNumberOfAgents(pSource) < 2)
    {
        double TregProbValue = 0.0;
        bool TregBool = TregThresholdRecNew(grid, pSource, TregProbValue);

        if (grid.getNumberOfAgents(pSource) == 1)
        {
            if (grid.hasAgentType(MAC, pSource))
                TregBool = g_Rand.getReal() < _PARAM(PARAM_TCELL_PROB_MOVE_TO_MAC);
            else
                TregBool = g_Rand.getReal() < _PARAM(PARAM_TCELL_PROB_MOVE_TO_TCELL);
        }

//        if (pSource == Pos (53,48))
//            std::cout << "Treg: " << pSource  << "  "  << TregProbValue << "  " << TregBool << std::endl;

        return (TregBool && g_Rand.getReal() < (_PARAM(PARAM_TREG_MAX_RECRUITMENT) * TregProbValue));
    }

    else
    {
        return false;
    }
}


