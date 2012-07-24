/*
 * grviz.h
 *
 *  Created on: Jan 20, 2010
 *      Author: mohammed
 */

#ifndef GRVIZ_H_
#define GRVIZ_H_

#include "simulation/gr.h"
#include <algorithm>
#include <vector>
#include <list>
#include <iostream>
#include <assert.h>
#include <math.h>
#include <float.h>

#ifdef _MSC_VER
	// with MSVC it is required to include <windows.h> prior to including OpenGL headers
	#define WIN32_LEAN_AND_MEAN
	#include <windows.h>

	// now get rid of the ugly min and max macros that microsoft wants to force upon us,
	// as they conflict with std::min and std::max
	#undef min
	#undef max

	// get rid of warning C4351
	#pragma warning(disable : 4351)

	// get rid of silly warnings about sscanf being unsafe
	#define _CRT_SECURE_NO_WARNINGS

	// microsoft's math.h is not C99 compliant, round and roundf are missing
    // Fixed macro issue.
    template<typename T>
    inline T roundf(const T& x) { return (x - floor(x)) > 0.5f ? ceil(x) : floor(x); }
    template<typename T>
    inline T round(const T& x) { return (x - floor(x)) > 0.5f ? ceil(x) : floor(x); }
#endif

#ifdef __APPLE__
	#include <OpenGL/gl.h>
	#include <OpenGL/glu.h>
#else
	#include <GL/gl.h>
	#include <GL/glu.h>
#endif

#define _MIN_CLAMP_VALUE 0.0f
#define _MAX_CLAMP_VALUE 1.0f
#define _MIN_X 0.0f
#define _MAX_X 99.0f
#define _MIN_Y 0.0f
#define _MAX_Y 99.0f

/* GR border */
#define _AREA_THRESHOLD 0.5f
#define _GR_BORDER_SIZE 5.0f
#define _GR_BORDER_ALPHA 0.9f
#define _GR_BORDER_COLOR Qt::red
#define _SLIDER_GRANULOMA_BORDER_THRESHOLD 20.0f

/* Color map */
#define _SLIDER_COLORMAP_HUE 10.0f
#define _SLIDER_COLORMAP_SAT 10.0f
#define _SLIDER_COLORMAP_VAL 10.0f
#define _SLIDER_COLORMAP_ALPHA 20.0f
#define _COLORMAP_LEGEND_PRECISION 5
#define _COLORMAP_ALPHA 0.2f
#define _COLORMAP_HUE_DELTA 0.0f
#define _COLORMAP_SAT_DELTA 0.0f
#define _COLORMAP_VAL_DELTA 0.0f
#define _COLORMAP_NR_BANDS 512

/* Isolines */
#define _SLIDER_ISOLINES_MIN 20.0f
#define _SLIDER_ISOLINES_MAX 20.0f
#define _SLIDER_ISOLINES_LABEL_PRECISION 4
#define _ISOLINES_TARGET_VALUE 0.5f
#define _ISOLINES_LINE_WIDTH 1.0f

/* Height plot */
#define _SLIDER_HEIGHTPLOT_ALPHA 20.0f
#define _SLIDER_HEIGHTPLOT_LABEL_PRECISION_SCALING 7
#define _SLIDER_HEIGHTPLOT_LABEL_PRECISION_CLAMPING 3
#define _HEIGHTPLOT_MAX_HEIGHT 8.0f
#define _HEIGHTPLOT_GRID_HEIGHT -4.0f
#define _HEIGHTPLOT_GRID_ALPHA 0.1f

/* Outcome */
#define _DEFAULT_OUTCOME_METHOD OUTCOME_NONE
#define _OUTCOME_SAMPLE_PERIOD 1
#define _OUTCOME_TEST_PERIOD TIME_STEPS_PER_DAY
#define _OUTCOME_ALPHA 0.05

/* Output */
#define _SNAPSHOT_PIC_PERIOD TIME_STEPS_PER_DAY
#define _SNAPSHOT_CSV_PERIOD 1

/* Simulation */
#define _DAYS_TO_SIMULATE 200
#define _TIMESTEPS_TO_SIMULATE (_DAYS_TO_SIMULATE * TIME_STEPS_PER_DAY)

/* Visualization */
#define _DRAW_AGENTS true
#define _DRAW_SMOKE true
#define _DRAW_GLYPHS false
#define _DRAW_ISOLINES false
#define _DRAW_HEIGHTPLOT false
#define _DRAW_GRANULOMA_BORDER true
#define _PRINT_TIME true
#define _PRINT_OUTCOME true

#endif /* GRVIZ_H_ */
