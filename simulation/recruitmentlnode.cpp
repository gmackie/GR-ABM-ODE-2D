/*
 * recruitmentlnode.cpp
 *
 *  Created on: Jun 17, 2010
 *      Author: mohammed
 */

#include "recruitmentlnode.h"
#include <stdlib.h>
#include <fstream>

RecruitmentLnODE::RecruitmentLnODE(const std::string& odeApp, const std::string& odeTmpFile)
	: _odeApp(odeApp)
	, _odeTmpFile(odeTmpFile)
{
	init();
}

RecruitmentLnODE::~RecruitmentLnODE()
{
}

void RecruitmentLnODE::init()
{
	_tcellQueue.clear();

	_tcellQueueCount[TCELL_TYPE_REG] = 0;
	_tcellQueueCount[TCELL_TYPE_GAM] = 0;
	_tcellQueueCount[TCELL_TYPE_CYT] = 0;

	_tcellTable[TCELL_TYPE_REG] = 0;
	_tcellTable[TCELL_TYPE_GAM] = 0;
	_tcellTable[TCELL_TYPE_CYT] = 0;

	_prevMiMci = 0;

	for (int i = 0; i < _nrConditions; i++)
	{
		switch (i)
		{
		case _idxNaiveCD4: // Naive CD4, used to be 1e5
			_odeInitialConditions[i] = _PARAM(PARAM_initN4);
			break;
		case _idxNaiveCD8: // Naive CD8, used to be .84e5
			_odeInitialConditions[i] = _PARAM(PARAM_initN8);
			break;
		default:
			_odeInitialConditions[i] = 0;
		}
	}
}

void RecruitmentLnODE::solveODE(GrStat&)
{
	std::string cmd = _odeApp;
	cmd += " '[";

	for (int i = 0; i < _nrConditions; i++)
	{
		char buf[1024];
		sprintf(buf, "%f,", _odeInitialConditions[i]);

		cmd += buf;
	}

	cmd += "]'";
	cmd += " > " + _odeTmpFile;

	assert_res(system(cmd.c_str()) == 0);

	std::ifstream inFile(_odeTmpFile.c_str());

	for (int i = 0; i < _nrConditions; i++)
	{
		assert(!inFile.eof());

		inFile >> _odeInitialConditions[i];
	}

	inFile.close();
}

void RecruitmentLnODE::updateQueue(GrStat& stats)
{
	// Compute the fluxes
	double tcellFlux[TCELL_TYPE_COUNT];

	tcellFlux[TCELL_TYPE_GAM] = _PARAM(PARAM_scaling_LN) *
			(_odeInitialConditions[_idxEffectorT8] + _odeInitialConditions[_idxEffectorTH1]);

	tcellFlux[TCELL_TYPE_CYT] = _PARAM(PARAM_scaling_LN) * _odeInitialConditions[_idxCTL];

	tcellFlux[TCELL_TYPE_REG] = .1 * tcellFlux[TCELL_TYPE_GAM];

	// Update the T cell queue and update the table that contains the lower bounds
	for (int i = 0; i < TCELL_TYPE_COUNT; i++)
	{
		TcellType type = (TcellType) i;

		if (tcellFlux[type] > _tcellTable[type] + 1)
		{
			const int count = tcellFlux[type] -_tcellTable[type];
			for (int j = 0; j < count; j++)
			{
				_tcellQueue.push_back(type);
				_tcellQueueCount[type]++;
			}

			_tcellTable[type] += count;
		}
		else if (tcellFlux[type] < _tcellTable[type])
		{
			_tcellTable[type] = tcellFlux[type];
		}
	}

	// Update the statistics
	stats.setFluxTgam(tcellFlux[TCELL_TYPE_GAM]);
	stats.setFluxTcyt(tcellFlux[TCELL_TYPE_CYT]);
	stats.setFluxTreg(tcellFlux[TCELL_TYPE_REG]);

	/*for (int i = 0; i < TCELL_TYPE_COUNT; i++)
	{
		TcellType type = (TcellType) i;

		switch (type)
		{
		case TCELL_TYPE_GAM:
			std::cout << "Tgam: ";
			break;
		case TCELL_TYPE_CYT:
			std::cout << "Tcyt: ";
			break;
		case TCELL_TYPE_REG:
			std::cout << "Treg: ";
			break;
		default:
			break;
		}

		std::cout << tcellFlux[type] << "\t" << _tcellDiscretizedTable[type] << "\t" << _tcellQueue[type] << std::endl;
	}*/
}

void RecruitmentLnODE::updateInitialConditions(GrStat& stats)
{
	int newMiMci = stats.getNrOfMacInfected() + stats.getNrOfMacCInfected();
	int delta = newMiMci - _prevMiMci;

	// update MDC, if the delta > 0
	if (delta > 0)
	{
		_odeInitialConditions[_idxMDC] += _PARAM(PARAM_scaling_LUNG) * _PARAM(PARAM_scaling_MDC) * delta;
	}

	// update _prevMiMci
	_prevMiMci = newMiMci;

	// update stats
	stats.setMDC(_odeInitialConditions[_idxMDC]);
}

void RecruitmentLnODE::recruitMacsGetTcellSources(GrSimulation& sim, GrStat& stats,
	ThresholdGridCellPtrList tcellSources[TCELL_TYPE_COUNT])
{
	GridCellPtrList& sources = sim.getGrid().getSources();

	// Recruit macs and make ordered lists of source pairs
	for (GridCellPtrList::iterator it = sources.begin();
		it != sources.end();
		it++)
	{
		double threshold;
		GridCell* pSource = *it;

		// if the source is caseated continue
		if (pSource->isCaseated())
			continue;

		if (MacRecruitmentThreshold(pSource, threshold))
		{
			stats.incNrSourcesMac();

			if (pSource->getNumberOfAgents() < 2 && !pSource->hasMac())
			{
				stats.incNrSourcesActiveMac();
			}
		}

		// macrophage recruitment
		if (pSource->getNumberOfAgents() < 2 && !pSource->hasMac())
			recruitMac(sim, pSource);

		if (TgamRecruitmentThreshold(pSource, threshold))
		{
			stats.incNrSourcesTgam();

			if (pSource->getNumberOfAgents() < 2)
			{
				stats.incNrSourcesActiveTgam();
				tcellSources[TCELL_TYPE_GAM].push_back(std::make_pair(threshold, pSource));
			}
		}

		if (TcytRecruitmentThreshold(pSource, threshold))
		{
			stats.incNrSourcesTcyt();

			if (pSource->getNumberOfAgents() < 2)
			{
				stats.incNrSourcesActiveTcyt();
				tcellSources[TCELL_TYPE_CYT].push_back(std::make_pair(threshold, pSource));
			}
		}

		if (TregRecruitmentThreshold(pSource, threshold))
		{
			stats.incNrSourcesTreg();

			if (pSource->getNumberOfAgents() < 2)
			{
				stats.incNrSourcesActiveTreg();
				tcellSources[TCELL_TYPE_REG].push_back(std::make_pair(threshold, pSource));
			}
		}
	}

	// sort the T cell sources
	for (int i = 0; i < TCELL_TYPE_COUNT; i++)
	{
		TcellType type = (TcellType) i;

		tcellSources[type].sort();
		tcellSources[type].reverse();
	}
}

void RecruitmentLnODE::recruitTcells(GrSimulation& sim, GrStat& stats,
		ThresholdGridCellPtrList tcellSources[TCELL_TYPE_COUNT])
{
	std::vector<TcellType> newTcellQueue;

	// pick a T cell from the queue
	while (_tcellQueue.size())
	{
		int idx = g_Rand.getInt((int) _tcellQueue.size());

		TcellType type = _tcellQueue[idx];
		_tcellQueue.erase(_tcellQueue.begin() + idx);

		// pick a source
		GridCell* pSource = NULL;
		for (ThresholdGridCellPtrList::iterator it = tcellSources[type].begin();
				it != tcellSources[type].end(); it++)
		{
			GridCell* pCurSource = it->second;
			assert (pSource->isSource() && !pSource->isCaseated());

			if (pCurSource->getNumberOfAgents() < 2)
			{
				// recruit
				pSource = pCurSource;
				break;
			}
		}

		if (!pSource)
		{
			newTcellQueue.push_back(type);
		}
		else
		{
			_tcellQueueCount[type]--;
			pSource->incNrRecruitments();

			// recruit it
			switch (type)
			{
			case TCELL_TYPE_CYT:
				sim.createTcyt(pSource->getRow(), pSource->getCol(),
					sim.getTime() - g_Rand.getInt(_PARAM(PARAM_TCELL_AGE), 1), TCYT_ACTIVE);
				break;
			case TCELL_TYPE_GAM:
				sim.createTgam(pSource->getRow(), pSource->getCol(),
					sim.getTime() - g_Rand.getInt(_PARAM(PARAM_TCELL_AGE), 1), TGAM_ACTIVE);
				break;
			case TCELL_TYPE_REG:
				sim.createTreg(pSource->getRow(), pSource->getCol(),
					sim.getTime() - g_Rand.getInt(_PARAM(PARAM_TCELL_AGE), 1), TREG_ACTIVE);
				break;
			default:
				break;
			}
		}
	}

	_tcellQueue = newTcellQueue;

	// update stats
	stats.setNrTgamQueued(_tcellQueueCount[TCELL_TYPE_GAM]);
	stats.setNrTcytQueued(_tcellQueueCount[TCELL_TYPE_CYT]);
	stats.setNrTregQueued(_tcellQueueCount[TCELL_TYPE_REG]);
}

void RecruitmentLnODE::recruit(GrSimulation& sim)
{
	/* Solve the ODE */
	GrStat& stats = sim.getStats();

	// update initial conditions (MDC)
	updateInitialConditions(stats);

	// solve the ODE for 10 minutes
	solveODE(stats);

	// update T cell queue according to new fluxes
	updateQueue(stats);

	/* Perform the actual recruitment */
	ThresholdGridCellPtrList tcellSources[TCELL_TYPE_COUNT];

	// recruit Macs
	recruitMacsGetTcellSources(sim, stats, tcellSources);

	// recruit T cells
	recruitTcells(sim, stats, tcellSources);
}

void RecruitmentLnODE::recruitMac(GrSimulation& sim, GridCell* pSource)
{
	assert(!pSource->isCaseated() && pSource->getNumberOfAgents() < 2 && !pSource->hasMac());

	// if the number of macrophages on the grid is less than _INITIAL_NUMBER_OF_MACROPHAGES,
	// recruit a resting macrophage
	if (sim.getStats().getNrOfMac() < _PARAM(PARAM_MAC_INIT_NUMBER))
	{
		sim.createMac(pSource->getRow(), pSource->getCol(),
			sim.getTime() - g_Rand.getInt(_PARAM(PARAM_MAC_AGE)), MAC_RESTING, false, false);
	}
	else
	{
		bool macThreshold = MacRecruitmentThreshold(pSource);

		if (macThreshold && g_Rand.getReal() < _PARAM(PARAM_MAC_PROB_RECRUITMENT))
		{
			pSource->incNrRecruitments();
			sim.createMac(pSource->getRow(), pSource->getCol(),
				sim.getTime() - g_Rand.getInt(_PARAM(PARAM_MAC_AGE)), MAC_RESTING, false, false);
		}
	}
}
