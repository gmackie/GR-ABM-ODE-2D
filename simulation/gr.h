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
#include "pos.h"
#include "params.h"
#include <valarray>
#include <boost/mpl/vector.hpp>

#ifndef SVN_VERSION
#define GR_VERSION "12022010"
#else
#define GR_VERSION "r" SVN_VERSION
#endif

//#ifndef __DIM__
//  #define __DIM__ (100)
//#endif //__DIM__
//#define NROWS (__DIM__)
//#define NCOLS (__DIM__)

//#define MOD_ROW(val) ((((val) + NROWS)) % NROWS)
//#define MOD_COL(val) ((((val) + NCOLS)) % NCOLS)
#define NOUTCOMES 2
#define TIME_STEPS_PER_DAY 144

// The number of compartments in a Moore neighborhood of a compartment.
#define MOORE_COUNT 9
#define MOORE_COUNT_DBL 9.0

#define PI 3.141592653

// Avogadro's Number
#define NAV 6.02e23

// used for conversion of conc. unit (M -> #/cell) based on cell and microcompartment volumes
#define DENSITY 1.25e11

// volume of a grid compartment in L
#define VOL 8.0e-12

//ratio of cytoplasmic to nuclear volume
#define KV 5.0

#define KDEG 4.58e-4

// molecular weight of IL10 in g/mol
#define MW_IL10 18600

// Number of seconds in a day
#define SECONDS_PER_DAY 86400

// Approximate limit for numerical accuracy solving the ODEs
#define MOLECULAR_ACCURACY 30.0

// Time step of agent movement in seconds
#define AGENT_TIME_STEP 600

// Number of Sig Figs of ODE Solver (Scaled to Power - i.e. 6 sig figs is 1e6)
#define ABS_TOL 1000000

#ifdef _DEBUG
#define assert_res(a) (assert(a))
#else
#define assert_res(a) (a)
#endif

// forward class declarations
class GrSimulation;
class GrGrid;
class Params;
class Tcell;
class Mac;
class Tgam;
class Tcyt;
class Treg;
class Agent;
class RecruitmentBase;
class Stats;

// typedefs
typedef std::list<Mac> MacList;
typedef std::list<Tgam> TgamList;
typedef std::list<Tcyt> TcytList;
typedef std::list<Treg> TregList;
typedef std::list<Mac*> MacPtrList;
typedef std::list<Tgam*> TgamPtrList;
typedef std::list<Tcyt*> TcytPtrList;
typedef std::list<Treg*> TregPtrList;
typedef std::vector<Mac> MacVector;
typedef std::vector<Tgam> TgamVector;
typedef std::vector<Tcyt> TcytVector;
typedef std::vector<Treg> TregVector;
typedef std::vector<Mac*> MacPtrVector;
typedef std::vector<Tgam*> TgamPtrVector;
typedef std::vector<Tcyt*> TcytPtrVector;
typedef std::vector<Treg*> TregPtrVector;
typedef std::pair<double, Pos> ThresholdPosPair;
typedef std::list<ThresholdPosPair> ThresholdPosList;

typedef enum {DIFF_REC_EQ = 0, DIFF_SOR_CORRECT = 1, DIFF_SOR_WRONG = 2, DIFF_REC_EQ_SWAP = 3, DIFF_ADE_SWAP = 4} DiffusionMethod;
typedef enum {OUTCOME_AREA = 0, OUTCOME_MTB = 1, OUTCOME_NONE = 2} OutcomeMethod;

typedef enum {MAC, TGAM, TCYT, TREG, NAGENTS} AgentType;
typedef boost::mpl::vector<Mac, Tgam, Tcyt, Treg> AgentTypes;

enum TcellType {TCELL_TYPE_CYT, TCELL_TYPE_REG, TCELL_TYPE_GAM, TCELL_TYPE_COUNT};
typedef int State;

unsigned int createTimeSeed();
// global variables
extern Rand g_Rand;

#endif /* GR_H_ */
