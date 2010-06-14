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

#define GR_VERSION "12022010"
#define NROWS 100
#define NCOLS 100
#define MOD_ROW(val) ((((val) + NROWS)) % NROWS)
#define MOD_COL(val) ((((val) + NCOLS)) % NCOLS)
#define NOUTCOMES 2

// The number of compartments in a Moore neighborhood of a compartment.
#define MOORE_COUNT 9

#ifdef _DEBUG
#define assert_res(a) (assert(a))
#else
#define assert_res(a) (a)
#endif

// forward class declarations
class GridCell;
class GrGrid;
class Params;
class Tcell;
class Mac;
class Tgam;
class Tcyt;
class Treg;
class Agent;

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

typedef enum {DIFF_REC_EQ = 0, DIFF_SOR_CORRECT = 1, DIFF_SOR_WRONG = 2} DiffusionMethod;
typedef enum {OUTCOME_AREA = 0, OUTCOME_MTB = 1, OUTCOME_NONE = 2} OutcomeMethod;

typedef enum {MAC, TGAM, TCYT, TREG} AgentType;
typedef enum {MAC_DEAD, MAC_RESTING, MAC_INFECTED, MAC_CINFECTED, MAC_ACTIVE} MacState;
typedef enum {TGAM_DEAD, TGAM_ACTIVE, TGAM_DOWN_REGULATED} TgamState;
typedef enum {TCYT_DEAD, TCYT_ACTIVE, TCYT_DOWN_REGULATED} TcytState;
typedef enum {TREG_DEAD, TREG_ACTIVE} TregState;


// global variables
extern Rand g_Rand;

#endif /* GR_H_ */
