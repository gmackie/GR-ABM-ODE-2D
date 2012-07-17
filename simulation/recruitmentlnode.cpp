/*
 * recruitmentlnode.cpp
 *
 *  Created on: Jun 17, 2010
 *      Author: mohammed
 */

#include "recruitmentlnode.h"
#include "serialization.h"
#include "grsimulation.h"
#include <stdlib.h>
#include <fstream>

const std::string RecruitmentLnODE::_ClassName = "RecruitmentLnODE";

RecruitmentLnODE::RecruitmentLnODE(const std::string& odeApp, const std::string& odeTmpFile, std::istream& in)
	: _odeApp(odeApp)
	, _odeTmpFile(odeTmpFile)
{
	RecruitmentLnODE::deserialize(in);
}

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

void RecruitmentLnODE::solveODE(const int, const Stats&, Stats&)
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

void RecruitmentLnODE::updateQueue(const int time, Stats& stats)
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
			const int count = (int) (tcellFlux[type] -_tcellTable[type]);
			for (int j = 0; j < count; j++)
			{
				_tcellQueue.push_back(TcellTypePair(time - g_Rand.getInt(_PARAM(PARAM_TCELL_AGE), 1), type));
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
	stats.setFlux<Tgam>(int(tcellFlux[TCELL_TYPE_GAM]));
	stats.setFlux<Tcyt>(int(tcellFlux[TCELL_TYPE_CYT]));
	stats.setFlux<Treg>(int(tcellFlux[TCELL_TYPE_REG]));

	// update stats
	stats.setNrQueued<Tgam>(_tcellQueueCount[TCELL_TYPE_GAM]);
	stats.setNrQueued<Tcyt>(_tcellQueueCount[TCELL_TYPE_CYT]);
	stats.setNrQueued<Treg>(_tcellQueueCount[TCELL_TYPE_REG]);

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
bool compare_tcell_pos(ThresholdPosPair first, ThresholdPosPair second) {
      return first.first < second.first;    //Keep deterministic by not sorting on pointers
}
void RecruitmentLnODE::updateInitialConditions(Stats& stats)
{
	int newMiMci = stats.getNrOfMacs(Mac::MAC_INFECTED) + stats.getNrOfMacs(Mac::MAC_CINFECTED);
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

void RecruitmentLnODE::recruitMacsGetTcellSources(GrSimulation& sim, Stats& stats,
	ThresholdPosList tcellSources[TCELL_TYPE_COUNT])
{
  GrGrid& grid = sim.getGrid();
	const PosVector& sources = sim.getGrid().getSources();

	// Recruit macs and make ordered lists of source pairs
	for (PosVector::const_iterator it = sources.begin();
		it != sources.end();
		it++)
	{
		double threshold;

		// if the source is caseated continue
		if (grid.isCaseated(*it))
			continue;

		if (MacRecruitmentThreshold(grid, *it, threshold))
		{
			++stats.getNrSources<Mac>();

			if (grid.getNumberOfAgents(*it) < 2 && !grid.hasAgentType(MAC, *it))
			{
				++stats.getNrSourcesActive<Mac>();
			}
			else
			{
				++stats.getNrSourcesCrowded<Mac>();
			}
		}

		// macrophage recruitment
		if (grid.getNumberOfAgents(*it) < 2 && !grid.hasAgentType(MAC, *it))
			recruitMac(sim, *it);

		if (TgamRecruitmentThreshold(grid, *it, threshold))
		{
			++stats.getNrSources<Tgam>();

			if (grid.getNumberOfAgents(*it) < 2)
			{
				++stats.getNrSourcesActive<Tgam>();
				tcellSources[TCELL_TYPE_GAM].push_back(std::make_pair(threshold, *it));
			}
			else
			{
				++stats.getNrSourcesCrowded<Tgam>();
			}
		}

		if (TcytRecruitmentThreshold(grid, *it, threshold))
		{
			++stats.getNrSources<Tcyt>();

			if (grid.getNumberOfAgents(*it) < 2)
			{
				++stats.getNrSourcesActive<Tcyt>();
				tcellSources[TCELL_TYPE_CYT].push_back(std::make_pair(threshold, *it));
			}
			else
			{
				++stats.getNrSourcesCrowded<Tcyt>();
			}
		}

		if (TregRecruitmentThreshold(grid, *it, threshold))
		{
			++stats.getNrSources<Treg>();

			if (grid.getNumberOfAgents(*it) < 2)
			{
				++stats.getNrSourcesActive<Treg>();
				tcellSources[TCELL_TYPE_REG].push_back(std::make_pair(threshold, *it));
			}
			else
			{
				++stats.getNrSourcesCrowded<Treg>();
			}
		}
	}

	// sort the T cell sources
	for (int i = 0; i < TCELL_TYPE_COUNT; i++)
	{
		TcellType type = (TcellType) i;

		tcellSources[type].sort(compare_tcell_pos);
		tcellSources[type].reverse();
	}
}

void RecruitmentLnODE::recruitTcells(GrSimulation& sim, Stats& stats,
		ThresholdPosList tcellSources[TCELL_TYPE_COUNT])
{
	std::vector<TcellTypePair> newTcellQueue;

	// pick a T cell from the queue
	while (_tcellQueue.size())
	{
		int idx = g_Rand.getInt((int) _tcellQueue.size());

		TcellTypePair tcell = _tcellQueue[idx];
		_tcellQueue.erase(_tcellQueue.begin() + idx);

		// check whether the T cell has died while in the queue
		if (tcell._birthtime + _PARAM(PARAM_TCELL_AGE) < sim.getTime())
		{
			++stats.getNrQueuedDie((AgentType)(tcell._type + MAC));
			continue;
		}

		// pick a source
		Pos* p = NULL;
		for (ThresholdPosList::iterator it = tcellSources[tcell._type].begin();
				it != tcellSources[tcell._type].end(); it++)
		{
			if (sim.getGrid().getNumberOfAgents(it->second) < 2)
			{
				// recruit
        p = &(it->second);
				break;
			}
		}

		if (!p)
		{
			newTcellQueue.push_back(tcell);
		}
		else
		{
			_tcellQueueCount[tcell._type]--;
			sim.getGrid().nRecruitments(*p)++;

			// recruit it
			switch (tcell._type)
			{
			case TCELL_TYPE_CYT:
				sim.createTcyt(p->x, p->y, tcell._birthtime, Tcyt::TCYT_ACTIVE);
				++stats.getNrRecruited<Tcyt>();
				break;
			case TCELL_TYPE_GAM:
				sim.createTgam(p->x, p->y, tcell._birthtime, Tgam::TGAM_ACTIVE);
				++stats.getNrRecruited<Tgam>();
				break;
			case TCELL_TYPE_REG:
				sim.createTreg(p->x, p->y, tcell._birthtime, Treg::TREG_ACTIVE);
				++stats.getNrRecruited<Treg>();
				break;
			default:
				assert(false);
				break;
			}
		}
	}


	_tcellQueue = newTcellQueue;
}

void RecruitmentLnODE::recruit(GrSimulation& sim)
{
	/* Solve the ODE */
	const Stats& statsPrevious = sim.getStatsPrevious();
	Stats& stats = sim.getStats();

	// update initial conditions (MDC)
	updateInitialConditions(stats);

	// solve the ODE for 10 minutes
	solveODE(sim.getTime(), statsPrevious, stats);

	// update T cell queue according to new fluxes
	updateQueue(sim.getTime(), stats);

	/* Perform the actual recruitment */
	ThresholdPosList tcellSources[TCELL_TYPE_COUNT];

	// recruit Macs
	recruitMacsGetTcellSources(sim, stats, tcellSources);

	// recruit T cells
	recruitTcells(sim, stats, tcellSources);
}

void RecruitmentLnODE::recruitMac(GrSimulation& sim, const Pos& p)
{
	assert(!sim.getGrid().isCaseated(p) && sim.getGrid().getNumberOfAgents(p) < 2 && !sim.getGrid().hasAgentType(MAC,p));

	// if the number of macrophages on the grid is less than _INITIAL_NUMBER_OF_MACROPHAGES,
	// recruit a resting macrophage
	if (sim.getStats().getNrOfMacs() < _PARAM(PARAM_MAC_INIT_NUMBER))
	{
		Mac* newMac = sim.createMac(p.x, p.y,
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
		bool macThreshold = MacRecruitmentThreshold(sim.getGrid(), p);

		if (macThreshold && g_Rand.getReal() < _PARAM(PARAM_MAC_PROB_RECRUITMENT))
		{
			sim.getGrid().nRecruitments(p)++;
  		Mac* newMac = sim.createMac(p.x, p.y,
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

void RecruitmentLnODE::serialize(std::ostream& out) const
{
	assert(out.good());

	Serialization::writeHeader(out, RecruitmentLnODE::_ClassName);

	for (int i = 0; i < TCELL_TYPE_COUNT; ++i)
	{
		out << _tcellTable[i] << std::endl;
	}

	for (int i = 0; i < TCELL_TYPE_COUNT; ++i)
	{
		out << _tcellQueueCount[i] << std::endl;
	}

	out << _tcellQueue.size()  << std::endl;
	for (std::vector<TcellTypePair>::const_iterator it = _tcellQueue.begin(); it != _tcellQueue.end(); ++it)
	{
		out << it->_birthtime << " " << it->_type << std::endl;
	}

	out << _prevMiMci << std::endl;

	for (int i = 0; i < _nrConditions; ++i)
	{
		out << _odeInitialConditions[i] << std::endl;
	}

	// _odeApp and _odeTmpFile are defined in the constructor, so they are not written here.
	// They are run dependent: a simulation state may be saved on one machine and copied
	// to another, where the _odeApp and _odeTmpFile might be in different locations.

	Serialization::writeFooter(out, RecruitmentLnODE::_ClassName);
}

void RecruitmentLnODE::deserialize(std::istream& in)
{
	assert(in.good());

	Serialization::readHeader(in, RecruitmentLnODE::_ClassName);

	for (int i = 0; i < TCELL_TYPE_COUNT; ++i)
	{
		in >> _tcellTable[i];
	}

	for (int i = 0; i < TCELL_TYPE_COUNT; ++i)
	{
		in >> _tcellQueueCount[i];
	}

	std::vector<TcellTypePair>::size_type queueSize;
	in >> queueSize;
	for (std::vector<TcellTypePair>::size_type i = 0; i < queueSize; ++i)
	{
		int birthtime;
		int temp;
		TcellType type;

		in >> birthtime;
		in >> temp;
		type = (TcellType) temp;

		_tcellQueue.push_back(TcellTypePair(birthtime, type));
	}

	in >> _prevMiMci;

	for (int i = 0; i < _nrConditions; ++i)
	{
		in >> _odeInitialConditions[i];
	}

	// _odeApp and _odeTmpFile are defined in the constructor, so they are not read here.
	// They are run dependent: a simulation state may be saved on one machine and copied
	// to another, where the _odeApp and _odeTmpFile might be in different locations.

	Serialization::readFooter(in, RecruitmentLnODE::_ClassName);
}
