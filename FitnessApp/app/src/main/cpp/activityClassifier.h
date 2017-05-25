//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: activityClassifier.h
//
// MATLAB Coder version            : 3.1
// C/C++ source code generated on  : 22-Apr-2017 17:39:26
//
#ifndef ACTIVITYCLASSIFIER_H
#define ACTIVITYCLASSIFIER_H

// Include Files
#include <cmath>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "activityClassifier_types.h"

// Variable Declarations

// Function Declarations
extern double activityClassifier(const double accX[50], const double accY[50],
  const double accZ[50], const double gyrX[50], const double gyrY[50], const
  double gyrZ[50]);
extern void activityClassifier_initialize();

#endif

//
// File trailer for activityClassifier.h
//
// [EOF]
//
