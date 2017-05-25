//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: resultant.h
//
// MATLAB Coder version            : 3.1
// C/C++ source code generated on  : 25-Mar-2017 15:03:24
//
#ifndef RESULTANT_H
#define RESULTANT_H

// Include Files
#include <cmath>
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "resultant_types.h"

// Function Declarations
extern double resultantAcc(const double X[3]);
extern double minDataMean(const double X[1000], const int len);
extern double resultantGyro(const double X[3]);
extern void resultant_initialize();

#endif

//
// File trailer for resultant.h
//
// [EOF]
//
