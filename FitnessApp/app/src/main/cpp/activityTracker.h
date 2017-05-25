//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: activityTracker.h
//
// MATLAB Coder version            : 3.1
// C/C++ source code generated on  : 08-Apr-2017 10:04:52
//
#ifndef ACTIVITYTRACKER_H
#define ACTIVITYTRACKER_H

// Include Files
#include <cmath>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "activityTracker_types.h"

// Variable Declarations

// Function Declarations
//extern void activityTracker(const double X[600], double Y[5]);

extern void activityTracker(double X[600], double Y[5]);
extern void activityTracker_initialize();
extern void activityTracker_terminate();

#endif

//
// File trailer for activityTracker.h
//
// [EOF]
//
