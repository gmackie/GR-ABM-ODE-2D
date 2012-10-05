/*
 * simulation.h
 *
 *  Created on: 3-sep-2008
 *      Author: S030858
 */

#ifndef SIMULATION_H_
#define SIMULATION_H_

#include "simulation/grsimulation.h"
#include <QThread>
#include <QMutex>
#include <QAtomicInt>

class Simulation : public QThread
{
  Q_OBJECT

private:
  mutable QMutex _mutex;
  mutable QMutex _modelMutex;
  int _time;
  GrSimulation* _gr;
  GrSimulation* _backbuffer;
  QAtomicInt _delay;
  bool _stopFlag;
  QAtomicInt _timeStepsToSimulate;
  QAtomicInt _mtbClearance;
  double _areaThreshold;
  OutcomeMethod _outcomeMethod;

  bool stopCondition();

public:
  void update();
  /* The following methods are thread-safe */
  Simulation(const Pos& dim);
  virtual ~Simulation();
  void run();
  void stop();
  void lock() const;
  void unlock() const;
  void modelLock() const;
  void modelUnlock() const;
  int getTime() const;
  int getTimeToSimulate() const;
  GrSimulation& getSim();
  const Pos& getSize() const
  {
    return getGrGrid().getRange();
  }
  double getAreaThreshold() const;
  DiffusionMethod getDiffusionMethod() const;
  OutcomeMethod getOutcomeMethod(int index) const;
  void getOutcomeParameters(int index, int& samplePeriod, int& testPeriod, double& alpha) const;
  void setDiffusionMethod(DiffusionMethod method);
  void setDaysToSimulate(int days);
  void setTimeToSimulate(int steps);
  void setMtbClearance(bool enable);
  void setDelay(int delay);
  void setAreaThreshold(double areaThreshold);
  void setOutcomeMethod(int index, OutcomeMethod method, double alpha, int testPeriod, int samplePeriod);

  /* The following methods are NOT thread-safe, first a lock must be obtained */
  const GrGrid& getGrGrid() const;
  const Stats& getStats() const;
  const MacList& getMacList() const;
  const TgamList& getTgamList() const;
  const TcytList& getTcytList() const;
  const TregList& getTregList() const;
  void loadState(std::istream& in);
  void saveState(std::ostream& out) const;
  void setRecruitmentMethod(RecruitmentMethod recrMethod);

  void setTnfrDynamics(bool tnfrDynamics);
  void setNfkbDynamics(bool nfkbDynamics);
  void setAdaptive(bool adaptive);
  void setTnfDepletionTimeStep(int tnfDepletionTimeStep);
  void setIl10DepletionTimeStep(int il10DepletionTimeStep);

  static QString getTimeStr(int simTime, int time);

public slots:
  void step();

signals:
  void stopConditionMet();
  void updated();
};

inline void Simulation::setOutcomeMethod(int index,
    OutcomeMethod method, double alpha, int testPeriod, int samplePeriod)
{
  lock();
  modelLock();
  _gr->setOutcomeMethod(index, method, alpha, testPeriod, samplePeriod);
  _backbuffer->setOutcomeMethod(index, method, alpha, testPeriod, samplePeriod);
  modelUnlock();
  unlock();
}

inline void Simulation::getOutcomeParameters(int index,
    int& samplePeriod, int& testPeriod, double& alpha) const
{
  lock();
  _backbuffer->getOutcomeParameters(index, samplePeriod, testPeriod, alpha);
  unlock();
}

inline double Simulation::getAreaThreshold() const
{
  double res;

  lock();
  res = _backbuffer->getAreaThreshold();
  unlock();

  return res;
}

inline void Simulation::setAreaThreshold(double areaThreshold)
{
  lock();
  modelLock();
  _gr->setAreaThreshold(areaThreshold);
  _backbuffer->setAreaThreshold(areaThreshold);
  modelUnlock();
  unlock();
}

inline void Simulation::setDelay(int delay)
{
  _delay = delay;
}

inline void Simulation::setDaysToSimulate(int days)
{
  setTimeToSimulate(TIME_STEPS_PER_DAY*days);
}

inline void Simulation::setTimeToSimulate(int steps)
{
  _timeStepsToSimulate = steps;
}

inline void Simulation::setMtbClearance(bool enable)
{
  _mtbClearance = enable;
}

inline void Simulation::setDiffusionMethod(DiffusionMethod method)
{
  lock();
  modelLock();
  _gr->setDiffusionMethod(method);
  _backbuffer->setDiffusionMethod(method);
  modelUnlock();
  unlock();
}

inline int Simulation::getTime() const
{
  return _time;
}

inline int Simulation::getTimeToSimulate() const
{
  return _timeStepsToSimulate;
}

inline DiffusionMethod Simulation::getDiffusionMethod() const
{
  DiffusionMethod res;

  lock();
  res = _backbuffer->getDiffusionMethod();
  unlock();

  return res;
}

inline OutcomeMethod Simulation::getOutcomeMethod(int index) const
{
  OutcomeMethod res;

  lock();
  res = _backbuffer->getOutcomeMethod(index);
  unlock();

  return res;
}

inline const Stats& Simulation::getStats() const
{
  return _backbuffer->getStats();
}

inline const MacList& Simulation::getMacList() const
{
  return _backbuffer->getMacList();
}

inline const TgamList& Simulation::getTgamList() const
{
  return _backbuffer->getTgamList();
}

inline const TcytList& Simulation::getTcytList() const
{
  return _backbuffer->getTcytList();
}

inline const TregList& Simulation::getTregList() const
{
  return _backbuffer->getTregList();
}

inline const GrGrid& Simulation::getGrGrid() const
{
  return _backbuffer->getGrid();
}

inline void Simulation::setRecruitmentMethod(RecruitmentMethod recrMethod)
{
  lock();
  modelLock();
  _gr->setRecruitmentMethod(recrMethod);
  _backbuffer->setRecruitmentMethod(recrMethod);
  unlock();
  modelUnlock();
}

inline void Simulation::setTnfrDynamics(bool tnfrDynamics)
{
  lock();
  modelLock();
  _gr->setTnfrDynamics(tnfrDynamics);
  _backbuffer->setTnfrDynamics(tnfrDynamics);
  unlock();
  modelUnlock();
}

inline void Simulation::setNfkbDynamics(bool nfkbDynamics)
{
  lock();
  modelLock();
  _gr->setNfkbDynamics(nfkbDynamics);
  _backbuffer->setNfkbDynamics(nfkbDynamics);
  unlock();
  modelUnlock();
}

inline void Simulation::setAdaptive(bool adaptive)
{
  lock();
  modelLock();
  _gr->setAdaptive(adaptive);
  _backbuffer->setAdaptive(adaptive);
  unlock();
  modelUnlock();
}

inline QString Simulation::getTimeStr(int simTime, int time)
{
  int days, hours, minutes;
  GrSimulation::convertSimTime(simTime, days, hours, minutes);

  QString str = QString("%1h%2m%3s - ").arg(time / (1000*60*60), 2, 10, QChar('0')).
                arg(((time % (60*60*1000)) / 60000), 2, 10, QChar('0')).
                arg((((time % (60*60*1000)) % 60000) / 1000), 2, 10, QChar('0')) +
                QString("%1d%2h%3m").arg(days, 3, 10, QChar('0')).
                arg(hours, 2, 10, QChar('0')).arg(minutes, 2, 10, QChar('0'));
  return str;
}

inline void Simulation::modelLock() const
{
  _modelMutex.lock();
}

inline void Simulation::modelUnlock() const
{
  _modelMutex.unlock();
}

inline GrSimulation& Simulation::getSim()
{
  return *_gr;
}

#endif /* SIMULATION_H_ */
