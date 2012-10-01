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
  moveToThread(this);
}

Simulation::~Simulation()
{
}

void Simulation::update()
{
  _updated = true;
  // get copies
  _grid = _gr->getGrid();
  _macList.clear();
  _time = _gr->getTime();
  const MacList& macList = _gr->getMacList();
  _macList.clear();
  for(typeof(macList.begin()) it = macList.begin(); it != macList.end(); it++)
    _macList.push_back(**it);
  const TgamList& tgamList = _gr->getTgamList();
  _tgamList.clear();
  for(typeof(tgamList.begin()) it = tgamList.begin(); it != tgamList.end(); it++)
    _tgamList.push_back(**it);
  const TcytList& tcytList = _gr->getTcytList();
  _tcytList.clear();
  for(typeof(tcytList.begin()) it = tcytList.begin(); it != tcytList.end(); it++)
    _tcytList.push_back(**it);
  const TregList& tregList = _gr->getTregList();
  _tregList.clear();
  for(typeof(tregList.begin()) it = tregList.begin(); it != tregList.end(); it++)
    _tregList.push_back(**it);
  _stats = _gr->getStats();
}

bool Simulation::stopCondition()
{
  return (_timeStepsToSimulate >= 0 && _time >= _timeStepsToSimulate) ||
         (_mtbClearance && _stats.getTotExtMtb() == 0 &&
          _stats.getTotIntMtb() == 0 &&
          _stats.getTotTNF() < DBL_EPSILON * 10.0);
}

void Simulation::step()
{
  _modelMutex.lock();
  _gr->solve();
  _modelMutex.unlock();

  lock();
  bool stop = stopCondition();

  if(stop)
    _gr->performT_Test();

  update();

  if(stop)
      emit stopConditionMet();

  unlock();
  msleep(_delay);

  if(!(stop || _stopFlag))
    QMetaObject::invokeMethod(this, "step", Qt::QueuedConnection);

}

void Simulation::run()
{
  QMetaObject::invokeMethod(this, "step", Qt::QueuedConnection);
  exec();
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
  _gr->load(in);
  update();
  _modelMutex.unlock();
}

void Simulation::saveState(std::ostream& out) const
{
  _modelMutex.lock();
  _gr->save(out);
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
