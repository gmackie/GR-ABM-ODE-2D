/*
 * gr.h
 *
 *  Created on: 02-nov-2009
 *      Author: M. El-Kebir
 */

#ifndef GR_H_
#define GR_H_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <list>
#include <vector>
#include <algorithm>
#include <utility>
#include "rand.h"

#ifndef SVN_VERSION
#define GR_VERSION "12022010"
#else
#define GR_VERSION "r" SVN_VERSION
#endif

#if defined(USE_FLOAT)
typedef float FLOAT_TYPE;
#define FLOAT_TYPE_PRECISION (FLT_DIG+2)
#else
typedef double FLOAT_TYPE;
#define FLOAT_TYPE_PRECISION (DBL_DIG+2)
#endif

#ifndef __DIM__
  #define __DIM__ (100)
#endif //__DIM__
#define NROWS (__DIM__)
#define NCOLS (__DIM__)

#define MOD_ROW(val) ((((val) + NROWS)) % NROWS)
#define MOD_COL(val) ((((val) + NCOLS)) % NCOLS)
#define NOUTCOMES 2
#define TIME_STEPS_PER_DAY 144

// The number of compartments in a Moore neighborhood of a compartment.
#define MOORE_COUNT 9
#define MOORE_COUNT_DBL 9.0

#define PI 3.141592653

#ifdef _DEBUG
#define assert_res(a) (assert(a))
#else
#define assert_res(a) (a)
#endif

// forward class declarations
class GrSimulation;
class GridCell;
class GrGrid;
class Params;
class Tcell;
class Mac;
class Tgam;
class Tcyt;
class Treg;
class Agent;
class RecruitmentBase;
class GrStats;

// typedefs
typedef std::list<GridCell*> GridCellPtrList;
typedef std::list<Mac> MacList;
typedef std::list<Tgam> TgamList;
typedef std::list<Tcyt> TcytList;
typedef std::list<Treg> TregList;
typedef std::list<Mac*> MacPtrList;
typedef std::list<Tgam*> TgamPtrList;
typedef std::list<Tcyt*> TcytPtrList;
typedef std::list<Treg*> TregPtrList;
typedef std::vector<GridCell*> GridCellPtrVector;
typedef std::vector<Mac> MacVector;
typedef std::vector<Tgam> TgamVector;
typedef std::vector<Tcyt> TcytVector;
typedef std::vector<Treg> TregVector;
typedef std::vector<Mac*> MacPtrVector;
typedef std::vector<Tgam*> TgamPtrVector;
typedef std::vector<Tcyt*> TcytPtrVector;
typedef std::vector<Treg*> TregPtrVector;
typedef std::pair<int, int> Pos;
typedef std::vector<Pos> PosVector;
typedef std::pair<double, GridCell*> ThresholdGridCellPtrPair;
typedef std::list<ThresholdGridCellPtrPair> ThresholdGridCellPtrList;

typedef enum {DIFF_REC_EQ = 0, DIFF_SOR_CORRECT = 1, DIFF_SOR_WRONG = 2, DIFF_REC_EQ_SWAP = 3} DiffusionMethod;
typedef enum {OUTCOME_AREA = 0, OUTCOME_MTB = 1, OUTCOME_NONE = 2} OutcomeMethod;
//typedef enum {RECR_PROB = 0, RECR_LN_ODE = 1} RecruitmentMethod;

typedef enum {MAC, TGAM, TCYT, TREG, NAGENTS} AgentType;
enum TcellType {TCELL_TYPE_CYT, TCELL_TYPE_REG, TCELL_TYPE_GAM, TCELL_TYPE_COUNT};
typedef enum {MAC_DEAD, MAC_RESTING, MAC_INFECTED, MAC_CINFECTED, MAC_ACTIVE, NMAC_STATES} MacState;
typedef enum {TGAM_DEAD, TGAM_ACTIVE, TGAM_DOWN_REGULATED,TGAM_ACTIVE_DOUBLE, TGAM_INDUCED_REG, NTGAM_STATES} TgamState;
typedef enum {TCYT_DEAD, TCYT_ACTIVE, TCYT_DOWN_REGULATED, NTCYT_STATES} TcytState;
typedef enum {TREG_DEAD, TREG_ACTIVE, NTREG_STATES} TregState;
typedef int State;

inline std::ostream& operator<<(std::ostream& s, const Pos& p) {
  return s<<'('<<p.first<<','<<p.second<<')';
}
inline std::istream& operator>>(std::istream& s, Pos& p) {
  char tmp;
  s>>tmp; assert(tmp == '(');
  s>>p.first;
  s>>tmp; assert(tmp == ',');
  s>>p.first;
  s>>tmp; assert(tmp == ')');
  return s;
}

// global variables
extern Rand g_Rand;

#endif /* GR_H_ */
