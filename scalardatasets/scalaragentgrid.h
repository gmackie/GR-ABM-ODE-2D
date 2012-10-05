/*
 * scalaragentgrid.h
 *
 *  Created on: 25-nov-2009
 *      Author: S030858
 */

#ifndef SCALARAGENTGRID_H_
#define SCALARAGENTGRID_H_

#include "simulation.h"
#include "scalaragentgridbase.h"

#define GET_BIT(x, i) (((x) >> (i)) & 1)
#define SET_BIT(i) ((1) << (i))

struct ScalarAgentItem
{
  int _bitMask;
  const Agent* _pAgent[2];
  int _nKillings;
  int _nRecruitments;
  int _nRecruitmentsMac;
  int _nRecruitmentsTgam;
  int _nRecruitmentsTcyt;
  int _nRecruitmentsTreg;
  int _nSecretions;
  double _attractant;
  double _TNF;
  double _CCL2;
  double _CCL5;
  double _CXCL9;
  double _shedTNFR2;
  double _il10;
  double _extMtb;
};

/**
 * reg dead cinf inf act res Treg Tcyt Tgam Mac cas src
 * 11  10   9    8   7   6   5    4    3    2   1   0
 */

class ScalarAgentGrid : public ScalarAgentGridBase
{
private:
  std::vector<ScalarAgentItem> _grid;
  std::vector<Mac> _macList;
  std::vector<Tgam> _tgamList;
  std::vector<Tcyt> _tcytList;
  std::vector<Treg> _tregList;
public:
  ScalarAgentGrid(size_t _DIM);
  ~ScalarAgentGrid();
  void evaluate(const Simulation* pSimulation);
  const std::vector<ScalarAgentItem>& getGrid() const;

  static const int _bitSrc = 0;
  static const int _bitCas = 1;
  static const int _bitMac = 2;
  static const int _bitTgam = 3;
  static const int _bitTcyt = 4;
  static const int _bitTreg = 5;
  static const int _bitMacResting = 6;
  static const int _bitMacActive = 7;
  static const int _bitMacInfected = 8;
  static const int _bitMacCInfected = 9;
  static const int _bitMacDead = 10;
  static const int _bitTgamActive = 11;
  static const int _bitTgamDownRegulated = 12;
  static const int _bitTgamDead = 13;
  static const int _bitTcytActive = 14;
  static const int _bitTcytDownRegulated = 15;
  static const int _bitTcytDead = 16;
  static const int _bitTregResting = 17;
  static const int _bitTregDead = 18;
  static const int _bitSrcMac = 19;
  static const int _bitSrcTgam = 20;
  static const int _bitSrcTcyt = 21;
  static const int _bitSrcTreg = 22;
  static const int _bitMacNFkB = 23;
  static const int _bitMacStat1 = 24;
  static const int _bitExtMtb = 25;
};

inline const std::vector<ScalarAgentItem>& ScalarAgentGrid::getGrid() const
{
  return _grid;
}

#endif /* SCALARAGENTGRID_H_ */
