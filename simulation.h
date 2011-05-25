/*
 * simulation.h
 *
 *  Created on: 3-sep-2008
 *      Author: S030858
 */

#ifndef SMOKESIMULATION_H_
#define SMOKESIMULATION_H_

#include "simulation/grsimulation.h"
#include <QThread>
#include <QMutex>

class Simulation : public QThread
{
	Q_OBJECT

private:
	mutable QMutex _mutex;
	mutable QMutex _modelMutex;
	int _time;
	GrSimulation _gr;
	GrGrid _grid;
	int _delay;
	bool _updated;
	bool _stopFlag;
	MacList _macList;
	TgamList _tgamList;
	TcytList _tcytList;
	TregList _tregList;
	GrStat _stats;
	int _timeStepsToSimulate;
	bool _mtbClearance;
	double _areaThreshold;
	OutcomeMethod _outcomeMethod;
	int _outcomeSamplePeriod;
	int _outcomeTestPeriod;
	double _outcomeAlpha;

	void update();
	bool stopCondition();

public:
	static const int _DIM = NROWS;

	/* The following methods are thread-safe */
	Simulation();
	virtual ~Simulation();
	void run();
	void stop();
	void lock() const;
	void unlock() const;
	void modelLock() const;
	void modelUnlock() const;
	int getTime() const;
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
	const GrStat& getStats() const;
	const MacList& getMacList() const;
	const TgamList& getTgamList() const;
	const TcytList& getTcytList() const;
	const TregList& getTregList() const;
	void loadState(std::istream& in);
	void saveState(std::ostream& out) const;
	void setRecruitment(RecruitmentBase* pRecruitment);
	
	void setTnfrDynamics(bool tnfrDynamics);

	static QString getTimeStr(int simTime, int time);

signals:
	void stopConditionMet();
};

inline void Simulation::setOutcomeMethod(int index,
	OutcomeMethod method, double alpha, int testPeriod, int samplePeriod)
{
	lock();
	_gr.setOutcomeMethod(index, method, alpha, testPeriod, samplePeriod);
	unlock();
}

inline void Simulation::getOutcomeParameters(int index,
	int& samplePeriod, int& testPeriod, double& alpha) const
{
	lock();
	_gr.getOutcomeParameters(index, samplePeriod, testPeriod, alpha);
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
	_gr.setAreaThreshold(areaThreshold);
	unlock();
}

inline void Simulation::setDelay(int delay)
{
	lock();
	_delay = delay;
	unlock();
}

inline void Simulation::setDaysToSimulate(int days)
{
	setTimeToSimulate(TIME_STEPS_PER_DAY*days);
}

inline void Simulation::setTimeToSimulate(int steps)
{
	lock();
	_timeStepsToSimulate = steps;
	unlock();
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
	_gr.setDiffusionMethod(method);
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
	int res;

	lock();
	res = _time;
	unlock();

	return res;
}

inline DiffusionMethod Simulation::getDiffusionMethod() const
{
	DiffusionMethod res;

	lock();
	res = _gr.getDiffusionMethod();
	unlock();

	return res;
}

inline OutcomeMethod Simulation::getOutcomeMethod(int index) const
{
	OutcomeMethod res;

	lock();
	res = _gr.getOutcomeMethod(index);
	unlock();

	return res;
}

inline const GrStat& Simulation::getStats() const
{
	return _stats;
}

inline const MacList& Simulation::getMacList() const
{
	return _macList;
}

inline const TgamList& Simulation::getTgamList() const
{
	return _tgamList;
}

inline const TcytList& Simulation::getTcytList() const
{
	return _tcytList;
}

inline const TregList& Simulation::getTregList() const
{
	return _tregList;
}

inline const GrGrid& Simulation::getGrGrid() const
{
	return _grid;
}

inline void Simulation::setRecruitment(RecruitmentBase* pRecruitment)
{
	_gr.setRecruitment(pRecruitment);
}

inline void Simulation::setTnfrDynamics(bool tnfrDynamics)
{
	_gr.setTnfrDynamics(tnfrDynamics);
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

#endif /* SMOKESIMULATION_H_ */
