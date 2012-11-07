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

/**
 * @brief Simulation thread that runs the model in the background.  Designed
 * to stop so the gui is able to render and account for *every* frame in order
 * to take accurate pictures and snapshots.
 */
class Simulation : public QThread
{
  Q_OBJECT

private:
  /// Locks the thread from updating the backbuffer
  mutable QMutex _mutex;
  /// Locks the thread from advancing the simulation.
  /// Useful if running settings need to change.
  mutable QMutex _modelMutex;
  /// Did the gui thread pick up the last simulation updated?
  /// If not, _updated is false and sim thread spins (not efficient)
  QAtomicInt _updated;
  /// Current time of simulation.
  int _time;
  /// Backup random number generator in dealing with cloning simulation.
  Rand rng;
  /// Running simulation
  GrSimulation* _gr;
  /// Updated simulation, gui thread uses this while holding _mutex
  GrSimulation* _backbuffer;
  /// #ms to wait before continuing to step
  QAtomicInt _delay;
  /// Been stopped (via condition)
  bool _stopFlag;
  /// Number of steps til stop
  QAtomicInt _timeStepsToSimulate;
  /// Stop on mtbClearance?
  QAtomicInt _mtbClearance;

  bool stopCondition();

public:
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
  bool getUpdated() const;
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
  void setUpdated(bool value=false);

  static QString getTimeStr(int simTime, int time);

public slots:
  void step();
  void update();

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

inline bool Simulation::getUpdated() const
{
  return _updated;
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

inline void Simulation::setUpdated(bool value)
{
  _updated = value;
}

#endif /* SIMULATION_H_ */
