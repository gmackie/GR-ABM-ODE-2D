/*
 * simulation.cpp
 *
 *  Created on: 3-sep-2008
 *      Author: S030858
 */

#include "simulation.h"
#include "gui/mainwindow.h"
#include <math.h>
#include <stdlib.h>
#include <cstring>

Simulation::Simulation(const Pos& dim)
	: _mutex(QMutex::Recursive)
	, _modelMutex(QMutex::Recursive)
	, _time(0)
	, _gr(new GrSimulation(dim))
	, _grid(dim)
	, _delay(0)
	, _updated(false)
	, _stopFlag(false)
	, _macList()
	, _tgamList()
	, _tcytList()
	, _tregList()
	, _stats()
	, _timeStepsToSimulate(_TIMESTEPS_TO_SIMULATE) //_DAYS_TO_SIMULATE _TIMESTEPS_TO_SIMULATE
	, _mtbClearance(true)
{
	_gr->init(false); // Molecular tracking not available in gui version of the model.
	_gr->setAreaThreshold(_AREA_THRESHOLD);
	update();
}

Simulation::~Simulation()
{
}

void Simulation::update()
{
	_updated = true;
	// get copies
	_grid = _gr->getGrid();
	_time = _gr->getTime();
	_macList = _gr->getMacList();
	_tgamList = _gr->getTgamList();
	_tcytList = _gr->getTcytList();
	_tregList = _gr->getTregList();
	_stats = _gr->getStats();
}

bool Simulation::stopCondition()
{
	return (_timeStepsToSimulate >= 0 && _time >= _timeStepsToSimulate) ||
			(_mtbClearance && _stats.getTotExtMtb() == 0 &&
			_stats.getTotIntMtb() == 0 &&
			_stats.getTotTNF() < DBL_EPSILON * 10.0);
}

void Simulation::run()
{
	lock();
	bool stop = stopCondition();
	unlock();
	while (!stop)
	{
		_modelMutex.lock();
		_gr->solve();
		_modelMutex.unlock();

		lock();
		update();
		stop = _stopFlag || stopCondition();
		unlock();

		msleep(_delay);
	}

	lock();
	if (stopCondition())
	{
		// Perform one last t-test so any display of simulation results
		// have current t-test results, not old ones (possibly very old ones
		// if the testing period is long).
		_gr->performT_Test();
		update(); // to get the latest t-test results from _gr
		emit stopConditionMet();
	}
	unlock();
}

void Simulation::lock() const
{
	_mutex.lock();
}

void Simulation::unlock() const
{
	_mutex.unlock();
}

void Simulation::stop()
{
	lock();
	_stopFlag = true;
	unlock();
}

void Simulation::loadState(std::istream& in)
{
	_modelMutex.lock();
	_stopFlag = false;
	_gr->deserialize(in);
	update();
	_modelMutex.unlock();
}

void Simulation::saveState(std::ostream& out) const
{
	_modelMutex.lock();
	_gr->serialize(out);
	_modelMutex.unlock();
}

void Simulation::setTnfDepletionTimeStep(int tnfDepletionTimeStep)
{
	lock();
	_gr->setTnfDepletionTimeStep(tnfDepletionTimeStep);
	unlock();
}

void Simulation::setIl10DepletionTimeStep(int il10DepletionTimeStep)
{
	lock();
	_gr->setIl10DepletionTimeStep(il10DepletionTimeStep);
	unlock();
}
