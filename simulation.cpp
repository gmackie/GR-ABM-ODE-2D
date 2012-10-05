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
  , _backbuffer(new GrSimulation(dim))
  , _delay(0)
  , _stopFlag(false)
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
  _backbuffer = _gr->clone(_backbuffer);  //Deep copy without deallocating
  _time = _backbuffer->getTime();
  emit updated();
}

bool Simulation::stopCondition()
{
  return (_timeStepsToSimulate >= 0 && _time >= _timeStepsToSimulate) ||
         (_mtbClearance && _gr->getStats().getTotExtMtb() == 0 &&
          _gr->getStats().getTotIntMtb() == 0 &&
          _gr->getStats().getTotTNF() < DBL_EPSILON * 10.0);
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
  lock();
  _backbuffer->save(out);
  unlock();
}

void Simulation::setTnfDepletionTimeStep(int tnfDepletionTimeStep)
{
  lock();
  modelLock();
  _gr->setTnfDepletionTimeStep(tnfDepletionTimeStep);
  _backbuffer->setTnfDepletionTimeStep(tnfDepletionTimeStep);
  unlock();
  modelUnlock();
}

void Simulation::setIl10DepletionTimeStep(int il10DepletionTimeStep)
{
  lock();
  modelLock();
  _gr->setIl10DepletionTimeStep(il10DepletionTimeStep);
  _backbuffer->setIl10DepletionTimeStep(il10DepletionTimeStep);
  unlock();
  modelUnlock();
}
