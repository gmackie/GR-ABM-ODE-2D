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

void RecruitmentProb::recruit(GrSimulation &sim, int time)
{
    GrGrid& grid = sim.getGrid();
    const std::vector<Pos>& sources = grid.getSources();
    Stats& stats = sim.getStats();

    for (PosVector::const_iterator it = sources.begin(); it != sources.end(); it++)
    {

        if (grid.getNumberOfAgents(*it) == (int)GrGrid::MAX_AGENTS_PER_CELL || grid.isCaseated(*it))
        {
            continue;
        }

        for(size_t i=0; i<(GrGrid::MAX_AGENTS_PER_CELL - grid.getNumberOfAgents(*it)); i++)
        {
            recruit(sim, *it, i);
        }
    }

    // Update Rates
    stats.getRecRate<Mac>() = (double)stats.getNrRec(MAC)/time;
    stats.getRecRate<Tgam>() = (double)stats.getNrRec(TGAM)/time;
    stats.getRecRate<Tcyt>() = (double)stats.getNrRec(TCYT)/time;
    stats.getRecRate<Treg>() = (double)stats.getNrRec(TREG)/time;
}


void RecruitmentProb::recruit(GrSimulation &sim, const Pos& pSource, size_t cellRecNum)
{
    GrGrid& grid = sim.getGrid();
    Stats& stats = sim.getStats();

    // Recruitment comes down to a bit of dilemma. Do we have a vascular source 'traffic jam' or a
    // supply issue. If we want to assume that there is always a cell to recruit at a vascular source we should use the commented
    // code below that adjusts probabilities based on the cells available for recruitment. That code casues a change in flux of cells
    // due to the adjusting probabilities

    double runningProbability = 0.0;
    const double randomNumber = g_Rand.getReal();
//    const double evenProb = 0.25;
    const double evenProb = 0.50;
    const double evenProbTcell = (1.0/3.0);
    const bool pmac = PossibleRecruitMac(sim, pSource, cellRecNum)  && !grid.hasAgentType(MAC, pSource);
    const bool ptgam = PossibleRecruitTgam(sim, pSource, cellRecNum) && sim.getTCellRecruitmentBegun();
    const bool ptcyt = PossibleRecruitTcyt(sim, pSource, cellRecNum) && sim.getTCellRecruitmentBegun();
    const bool ptreg = PossibleRecruitTreg(sim, pSource, cellRecNum) && sim.getTCellRecruitmentBegun();
    const bool ptcell = ptgam || ptcyt || ptreg;
//    std::cout << "Possible Mac: " << pmac  << "  Possible Tgam: " << ptgam << "  Possible Tcyt: " << ptcyt << "  Possible Treg: " << ptreg << std::endl;

    // Update Recruitment Stats - Where cells are actively being recruited from
    stats.getNrSourcesActive<Mac>() += pmac;
    stats.getNrSourcesActive<Tgam>() += ptgam;
    stats.getNrSourcesActive<Tcyt>() += ptcyt;
    stats.getNrSourcesActive<Treg>() += ptreg;

//    const size_t countPossible = pmac + ptgam + ptcyt + ptreg;

    const size_t countPossible = pmac + ptcell;


    if (pmac)
    {
        runningProbability += (evenProb * pmac) + (evenProb / countPossible) * (!ptcell);

        if ((randomNumber < runningProbability))
        {
            recruitCellMac(sim, pSource);
            ++stats.getNrRec(MAC);
//            std::cout << "Recruit Mac and Break" << std::endl;
            return;
        }

    }

    if (ptcell)
    {
        runningProbability += (evenProb * ptcell) + (evenProb / countPossible) * (!pmac);

        if ((randomNumber <= runningProbability))
        {

            double runningProbabilityTcell = 0.0;
            const size_t countPossibleTcell = ptgam + ptcyt + ptreg;
            const double randomNumberTcell = g_Rand.getReal();

                if (ptgam)
                {
                    runningProbabilityTcell += (evenProbTcell * ptgam) + (evenProbTcell/countPossibleTcell) * (!ptcyt + !ptreg);

                    if ((randomNumberTcell < runningProbabilityTcell))
                    {
                        recruitCellTgam(sim, pSource);
                        ++stats.getNrRec(TGAM);
            //            std::cout << "Recruit Tgam and Break" << std::endl;
                        return;
                    }

                }


                if (ptcyt)
                {
                    runningProbabilityTcell += (evenProbTcell * ptcyt) + (evenProbTcell / countPossibleTcell) * (!ptgam + !ptreg);

                    if ((randomNumberTcell < runningProbabilityTcell))
                    {
                        recruitCellTcyt(sim, pSource);
                        ++stats.getNrRec(TCYT);
            //            std::cout << "Recruit Tcyt and Break" << std::endl;
                        return;
                    }

                }


                if (ptreg)
                {
                    runningProbabilityTcell += (evenProbTcell * ptreg) + (evenProbTcell / countPossibleTcell) * (!ptgam + !ptcyt);

                    if ((randomNumberTcell <= runningProbabilityTcell))
                    {
                        recruitCellTreg(sim, pSource);
                        ++stats.getNrRec(TREG);
            //            std::cout << "Recruit Treg and Break" << std::endl;
                        return;
                    }

                }

            throw std::runtime_error("Something has gone horribly wrong with Tcell Recruitment -- Please check recruitment algorithm");
            return;
        }

    }

//    if (pmac)
//    {
//        runningProbability += (evenProb * pmac) + (evenProb / countPossible) * (!ptgam +
//                            !ptcyt + !ptreg);

//        if ((randomNumber < runningProbability))
//        {
//            recruitCellMac(sim, pSource);
//            ++stats.getNrRec(MAC);
////            std::cout << "Recruit Mac and Break" << std::endl;
//            return;
//        }

//    }

//    if (ptgam)
//    {
//        runningProbability += (evenProb * ptgam) + (evenProb / countPossible) * (!pmac +
//                            !ptcyt + !ptreg);

//        if ((randomNumber < runningProbability))
//        {
//            recruitCellTgam(sim, pSource);
//            ++stats.getNrRec(TGAM);
////            std::cout << "Recruit Tgam and Break" << std::endl;
//            return;
//        }

//    }

//    if (ptcyt)
//    {
//        runningProbability += (evenProb * ptcyt) + (evenProb / countPossible) * (!ptgam +
//                            !pmac + !ptreg);

//        if ((randomNumber < runningProbability))
//        {
//            recruitCellTcyt(sim, pSource);
//            ++stats.getNrRec(TCYT);
////            std::cout << "Recruit Tcyt and Break" << std::endl;
//            return;
//        }

//    }

//    if (ptreg)
//    {
//        runningProbability += (evenProb * ptreg) + (evenProb / countPossible) * (!ptgam +
//                            !ptcyt + !pmac);

//        if ((randomNumber < runningProbability))
//        {
//            recruitCellTreg(sim, pSource);
//            ++stats.getNrRec(TREG);
////            std::cout << "Recruit Treg and Break" << std::endl;
//            return;
//        }

//    }

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
    ++grid.nRecruitmentsMac(pSource);
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
    ++grid.nRecruitmentsTgam(pSource);
    sim.createTgam(pSource.x, pSource.y,
        sim.getTime() - g_Rand.getInt(_PARAM(PARAM_TCELL_AGE), 1), Tgam::TGAM_ACTIVE);
}


inline void RecruitmentProb::recruitCellTcyt(GrSimulation &sim, const Pos &pSource)
{
    GrGrid& grid = sim.getGrid();

    assert(!grid.isCaseated(pSource) && grid.getNumberOfAgents(pSource) < 2);

    ++grid.nRecruitments(pSource);
    ++grid.nRecruitmentsTcyt(pSource);
    sim.createTcyt(pSource.x, pSource.y,
        sim.getTime() - g_Rand.getInt(_PARAM(PARAM_TCELL_AGE), 1), Tcyt::TCYT_ACTIVE);
}

inline void RecruitmentProb::recruitCellTreg(GrSimulation &sim, const Pos &pSource)
{
    GrGrid& grid = sim.getGrid();

    assert(!grid.isCaseated(pSource) && grid.getNumberOfAgents(pSource) < 2);

    ++grid.nRecruitments(pSource);
    ++grid.nRecruitmentsTreg(pSource);
    sim.createTreg(pSource.x, pSource.y,
        sim.getTime() - g_Rand.getInt(_PARAM(PARAM_TCELL_AGE), 1), Treg::TREG_ACTIVE);
}


bool RecruitmentProb::PossibleRecruitMac(GrSimulation& sim, const Pos& pSource, size_t cellRecNum)
{
    GrGrid& grid = sim.getGrid();
    Stats& stats = sim.getStats();

    double macProbValue = 0.0;
    size_t testNum = 0;
    bool macBool = MacThresholdRecNew(grid, pSource, macProbValue);

    if (macBool && cellRecNum == testNum)
        stats.getNrSources<Mac>() += macBool;

    if (!grid.isCaseated(pSource) && grid.getNumberOfAgents(pSource) < 2 && !grid.hasAgentType(MAC, pSource))
    {
        // if the number of macrophages on the grid is less than _INITIAL_NUMBER_OF_MACROPHAGES,
        // recruit a resting macrophage
        if (sim.getStats().getNrOfMacs() < _PARAM(PARAM_MAC_INIT_NUMBER))
            return true;

        return (macBool && g_Rand.getReal() < (_PARAM(PARAM_MAC_MAX_RECRUITMENT) * macProbValue));
    }

    else
    {
        return false;
    }
}

bool RecruitmentProb::PossibleRecruitTgam(GrSimulation& sim, const Pos& pSource, size_t cellRecNum)
{
    GrGrid& grid = sim.getGrid();
    Stats& stats = sim.getStats();

    double TgamProbValue = 0.0;
    size_t testNum = 0;
    bool TgamBool = TgamThresholdRecNew(grid, pSource, TgamProbValue);

    if (TgamBool && cellRecNum == testNum)
        stats.getNrSources<Tgam>() += TgamBool;

    if (!grid.isCaseated(pSource) && grid.getNumberOfAgents(pSource) < 2)
    {
        if (grid.getNumberOfAgents(pSource) == 1)
        {
            if (grid.hasAgentType(MAC, pSource))
                TgamBool = g_Rand.getReal() < _PARAM(PARAM_TCELL_PROB_MOVE_TO_MAC);
            else
                TgamBool = g_Rand.getReal() < _PARAM(PARAM_TCELL_PROB_MOVE_TO_TCELL);
        }

        return (TgamBool && g_Rand.getReal() < (_PARAM(PARAM_TGAM_MAX_RECRUITMENT) * TgamProbValue));

    }
    else
    {
        return false;
    }
}

bool RecruitmentProb::PossibleRecruitTcyt(GrSimulation& sim, const Pos& pSource, size_t cellRecNum)
{
    GrGrid& grid = sim.getGrid();
    Stats& stats = sim.getStats();

    double TcytProbValue = 0.0;
    size_t testNum = 0;
    bool TcytBool = TcytThresholdRecNew(grid, pSource, TcytProbValue);

    if (TcytBool && cellRecNum == testNum)
        stats.getNrSources<Tcyt>() += TcytBool;

    if (!grid.isCaseated(pSource) && grid.getNumberOfAgents(pSource) < 2)
    {
        if (grid.getNumberOfAgents(pSource) == 1)
        {
            if (grid.hasAgentType(MAC, pSource))
                TcytBool = g_Rand.getReal() < _PARAM(PARAM_TCELL_PROB_MOVE_TO_MAC);
            else
                TcytBool = g_Rand.getReal() < _PARAM(PARAM_TCELL_PROB_MOVE_TO_TCELL);
        }

        return (TcytBool && g_Rand.getReal() < (_PARAM(PARAM_TCYT_MAX_RECRUITMENT) * TcytProbValue));
    }
    else
    {
        return false;
    }
}

bool RecruitmentProb::PossibleRecruitTreg(GrSimulation& sim, const Pos& pSource, size_t cellRecNum)
{
    GrGrid& grid = sim.getGrid();
    Stats& stats = sim.getStats();

    double TregProbValue = 0.0;
    size_t testNum = 0;
    bool TregBool = TregThresholdRecNew(grid, pSource, TregProbValue);

    if (TregBool && cellRecNum == testNum)
        stats.getNrSources<Treg>() += TregBool;

    if (!grid.isCaseated(pSource) && grid.getNumberOfAgents(pSource) < 2)
    {
        if (grid.getNumberOfAgents(pSource) == 1)
        {
            if (grid.hasAgentType(MAC, pSource))
                TregBool = g_Rand.getReal() < _PARAM(PARAM_TCELL_PROB_MOVE_TO_MAC);
            else
                TregBool = g_Rand.getReal() < _PARAM(PARAM_TCELL_PROB_MOVE_TO_TCELL);
        }

        return (TregBool && g_Rand.getReal() < (_PARAM(PARAM_TREG_MAX_RECRUITMENT) * TregProbValue));
    }

    else
    {
        return false;
    }
}


