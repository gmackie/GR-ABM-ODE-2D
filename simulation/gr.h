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
#include "version.h"
#include "scalar.h"
#include <valarray>
#include <boost/mpl/vector.hpp>

#if 1
#define likely(x) __builtin_expect((x), 1)
#define unlikely(x) __builtin_expect((x), 0)
#else
#define unlikely(x) (x)
#define likely(x) (x)
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

// Define the Cash-Karp Butcher Table
#define CK_B_2_1 (1.0/5.0)
#define CK_B_3_1 (3.0/40.0)
#define CK_B_3_2 (9.0/40.0)
#define CK_B_4_1 (3.0/10.0)
#define CK_B_4_2 (-9.0/10.0)
#define CK_B_4_3 (6.0/5.0)
#define CK_B_5_1 (-11.0/54.0)
#define CK_B_5_2 (5.0/2.0)
#define CK_B_5_3 (-70.0/27.0)
#define CK_B_5_4 (35.0/27.0)
#define CK_B_6_1 (1631.0/55296.0)
#define CK_B_6_2 (175.0/512.0)
#define CK_B_6_3 (575.0/13824.0)
#define CK_B_6_4 (44275.0/110592.0)
#define CK_B_6_5 (253.0/4096.0)
#define CK_C_1 (37.0/378.0)
#define CK_C_3 (250.0/621.0)
#define CK_C_4 (125.0/594.0)
#define CK_C_6 (512.0/1771.0)
#define CK_DC_1 (CK_C_1-(2825.0/27648.0))
#define CK_DC_3 (CK_C_3-(18575.0/48384.0))
#define CK_DC_4 (CK_C_4-(13525.0/55296.0))
#define CK_DC_5 (-277.0/14336.0)
#define CK_DC_6 (CK_C_6-(0.25))

// Define Error Estimate Factors
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
typedef std::vector<Mac*> MacList;
typedef std::vector<Tgam*> TgamList;
typedef std::vector<Tcyt*> TcytList;
typedef std::vector<Treg*> TregList;
typedef std::vector<Mac*> MacPtrList;
typedef std::vector<Tgam*> TgamPtrList;
typedef std::vector<Tcyt*> TcytPtrList;
typedef std::vector<Treg*> TregPtrList;
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
typedef enum {RECR_PROB = 0, RECR_LN_ODE_PROXY = 1, RECR_LN_ODE_PURE = 2, RECR_LN_ODE = 3} RecruitmentMethod;


typedef enum {MAC, TGAM, TCYT, TREG, NAGENTS} AgentType;
typedef boost::mpl::vector<Mac, Tgam, Tcyt, Treg> AgentTypes;

enum TcellType {TCELL_TYPE_CYT, TCELL_TYPE_REG, TCELL_TYPE_GAM, TCELL_TYPE_COUNT};
typedef int State;

unsigned int createTimeSeed();
// global variables
extern Rand g_Rand;

//Print agent type
inline std::ostream& operator<<(std::ostream& os, AgentType at) {
  switch(at) {
  case MAC: os<<"Mac"; break;
  case TGAM: os<<"Tgam"; break;
  case TCYT: os<<"Tcyt"; break;
  case TREG: os<<"Treg"; break;
  default: break;
  }
  return os;
}

#endif /* GR_H_ */
