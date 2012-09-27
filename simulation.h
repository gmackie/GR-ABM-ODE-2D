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
  QAtomicInt _time;
  GrSimulation* _gr;
  GrGrid _grid;
  QAtomicInt _delay;
  bool _updated;
  bool _stopFlag;
  std::vector<Mac> _macList;
  std::vector<Tgam> _tgamList;
  std::vector<Tcyt> _tcytList;
  std::vector<Treg> _tregList;
  Stats _stats;
  QAtomicInt _timeStepsToSimulate;
  bool _mtbClearance;
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
  bool getUpdated() const;
  double getAreaThreshold() const;
  DiffusionMethod getDiffusionMethod() const;
  OutcomeMethod getOutcomeMethod(int index) const;
  void getOutcomeParameters(int index, int& samplePeriod, int& testPeriod, double& alpha) const;
  void setUpdated(bool flag);
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
  const std::vector<Mac>& getMacList() const;
  const std::vector<Tgam>& getTgamList() const;
  const std::vector<Tcyt>& getTcytList() const;
  const std::vector<Treg>& getTregList() const;
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
};

inline void Simulation::setOutcomeMethod(int index,
    OutcomeMethod method, double alpha, int testPeriod, int samplePeriod)
{
  lock();
  _gr->setOutcomeMethod(index, method, alpha, testPeriod, samplePeriod);
  unlock();
}

inline void Simulation::getOutcomeParameters(int index,
    int& samplePeriod, int& testPeriod, double& alpha) const
{
  lock();
  _gr->getOutcomeParameters(index, samplePeriod, testPeriod, alpha);
  unlock();
}

inline double Simulation::getAreaThreshold() const
{
  double res;

  lock();
  res = _areaThreshold;
  unlock();

  return res;
}

inline void Simulation::setAreaThreshold(double areaThreshold)
{
  lock();
  _gr->setAreaThreshold(areaThreshold);
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
  lock();
  _mtbClearance = enable;
  unlock();
}

inline void Simulation::setDiffusionMethod(DiffusionMethod method)
{
  lock();
  _gr->setDiffusionMethod(method);
  unlock();
}

inline bool Simulation::getUpdated() const
{
  bool res;

  lock();
  res = _updated;
  unlock();

  return res;
}

inline void Simulation::setUpdated(bool flag)
{
  lock();
  _updated = flag;
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
  res = _gr->getDiffusionMethod();
  unlock();

  return res;
}

inline OutcomeMethod Simulation::getOutcomeMethod(int index) const
{
  OutcomeMethod res;

  lock();
  res = _gr->getOutcomeMethod(index);
  unlock();

  return res;
}

inline const Stats& Simulation::getStats() const
{
  return _stats;
}

inline const std::vector<Mac>& Simulation::getMacList() const
{
  return _macList;
}

inline const std::vector<Tgam>& Simulation::getTgamList() const
{
  return _tgamList;
}

inline const std::vector<Tcyt>& Simulation::getTcytList() const
{
  return _tcytList;
}

inline const std::vector<Treg>& Simulation::getTregList() const
{
  return _tregList;
}

inline const GrGrid& Simulation::getGrGrid() const
{
  return _grid;
}

inline void Simulation::setRecruitmentMethod(RecruitmentMethod recrMethod)
{
  _gr->setRecruitmentMethod(recrMethod);
}

inline void Simulation::setTnfrDynamics(bool tnfrDynamics)
{
  _gr->setTnfrDynamics(tnfrDynamics);
}

inline void Simulation::setNfkbDynamics(bool nfkbDynamics)
{
  _gr->setNfkbDynamics(nfkbDynamics);
}

inline void Simulation::setAdaptive(bool adaptive)
{
  _gr->setAdaptive(adaptive);
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
