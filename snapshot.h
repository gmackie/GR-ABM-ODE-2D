/*
 * snapshot.h
 *
 *  Created on: Dec 3, 2009
 *      Author: mohammed
 */

#ifndef SNAPSHOT_H_
#define SNAPSHOT_H_

#include <fstream>
#include <QString>

class QImage;
class Stats;
class Simulation;

class Snapshot
{
private:
  std::ofstream _outFile;
  const QString _dirName;

public:
  Snapshot(const QString& dirName, const QString& fileName);
  ~Snapshot();
  bool isGood() const;
  void takePicture(const int time, const QImage& image, int slice=-1, const QString prefix=QString());
  void takePicture(const int time, const QImage& image, const QString prefix);
  void takeSnapshot(const int time, const Stats& stats);
  void takeStateSnapshot(const int time, const Simulation& sim);
};

inline bool Snapshot::isGood() const
{
  return _outFile.good();
}

#endif /* SNAPSHOT_H_ */
