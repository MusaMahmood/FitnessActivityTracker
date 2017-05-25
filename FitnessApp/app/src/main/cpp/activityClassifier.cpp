//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: activityClassifier.cpp
//
// MATLAB Coder version            : 3.1
// C/C++ source code generated on  : 22-Apr-2017 17:39:26
//

// Include Files
#include "rt_nonfinite.h"
#include "activityClassifier.h"

// Type Definitions
#ifndef struct_emxArray__common
#define struct_emxArray__common

struct emxArray__common
{
  void *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif                                 //struct_emxArray__common

#ifndef struct_emxArray_boolean_T
#define struct_emxArray_boolean_T

struct emxArray_boolean_T
{
  boolean_T *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif                                 //struct_emxArray_boolean_T

#ifndef struct_emxArray_int32_T
#define struct_emxArray_int32_T

struct emxArray_int32_T
{
  int *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif                                 //struct_emxArray_int32_T

#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T

struct emxArray_real_T
{
  double *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif                                 //struct_emxArray_real_T

// Variable Definitions

// Function Declarations
static void activityFeatureExtractionAcc(const double sX[50], emxArray_real_T *F);
static void activityFeatureExtractionGyro(const double sX[50], emxArray_real_T
  *F);
static void assignOutputs(const double y[50], const double x[50], const
  emxArray_real_T *iPk, emxArray_real_T *YpkOut, emxArray_real_T *XpkOut);
static void b_abs(const emxArray_real_T *x, emxArray_real_T *y);
static void b_diff(const emxArray_real_T *x, emxArray_real_T *y);
static void b_findLocalMaxima(const double yTemp[50], emxArray_real_T *iPk,
  emxArray_real_T *iInflect);
static void b_merge(int idx[245], double x[245], int offset, int np, int nq, int
                    iwork[245], double xwork[245]);
static void b_merge_block(int idx[245], double x[245], int offset, int n, int
  preSortLevel, int iwork[245], double xwork[245]);
static void b_sign(emxArray_real_T *x);
static void b_sort(emxArray_real_T *x, int dim, emxArray_int32_T *idx);
static double b_std(const double varargin_1[50]);
static void c_findPeaksSeparatedByMoreThanM(const emxArray_real_T *iPk,
  emxArray_real_T *idx);
static void c_sort(double x[245], int idx[245]);
static void combinePeaks(const emxArray_real_T *iPk, const emxArray_real_T *iInf,
  emxArray_real_T *iPkOut);
static void diff(const double x[50], double y[49]);
static void do_vectors(const emxArray_real_T *a, const emxArray_real_T *b,
  emxArray_real_T *c, emxArray_int32_T *ia, emxArray_int32_T *ib);
static void emxEnsureCapacity(emxArray__common *emxArray, int oldNumel, int
  elementSize);
static void emxFree_boolean_T(emxArray_boolean_T **pEmxArray);
static void emxFree_int32_T(emxArray_int32_T **pEmxArray);
static void emxFree_real_T(emxArray_real_T **pEmxArray);
static void emxInit_boolean_T(emxArray_boolean_T **pEmxArray, int numDimensions);
static void emxInit_int32_T(emxArray_int32_T **pEmxArray, int numDimensions);
static void emxInit_int32_T1(emxArray_int32_T **pEmxArray, int numDimensions);
static void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions);
static void emxInit_real_T1(emxArray_real_T **pEmxArray, int numDimensions);
static double eps(double x);
static void findLocalMaxima(const double yTemp[49], emxArray_real_T *iPk,
  emxArray_real_T *iInflect);
static void findpeaks(const double Yin[49], emxArray_real_T *Ypk,
                      emxArray_real_T *Xpk);
static void getAllPeaks(const double y[50], emxArray_real_T *iPk,
  emxArray_real_T *iInf, emxArray_real_T *iInflect);
static void keepAtMostNpPeaks(emxArray_real_T *idx);
static double knnclassify(const emxArray_real_T *testSamples);
static double mean(const emxArray_real_T *x);
static void merge(emxArray_int32_T *idx, emxArray_real_T *x, int offset, int np,
                  int nq, emxArray_int32_T *iwork, emxArray_real_T *xwork);
static void merge_block(emxArray_int32_T *idx, emxArray_real_T *x, int offset,
  int n, int preSortLevel, emxArray_int32_T *iwork, emxArray_real_T *xwork);
static void mrdivide(const emxArray_real_T *A, const emxArray_real_T *B,
                     emxArray_real_T *y);
static void removePeaksBelowMinPeakHeight(const double Y[50], emxArray_real_T
  *iPk, double Ph);
static void removePeaksBelowThreshold(const double Y[50], emxArray_real_T *iPk,
  double Th);
static void resultant(const double X[150], double Y[50]);
static double rt_hypotd_snf(double u0, double u1);
static double skip_to_last_equal_value(int *k, const emxArray_real_T *x);
static void sort(emxArray_real_T *x, emxArray_int32_T *idx);
static void sortIdx(emxArray_real_T *x, emxArray_int32_T *idx);
static double trapz(const double x[50]);
static void xgeqp3(emxArray_real_T *A, emxArray_real_T *tau, emxArray_int32_T
                   *jpvt);
static double xnrm2(int n, const emxArray_real_T *x, int ix0);
static void xscal(int n, double a, emxArray_real_T *x, int ix0);

// Function Definitions

//
// featureExtraction Summary of this function goes here if I ever feel like
// writing one up.
//  sX = samples X input
// Arguments    : const double sX[50]
//                emxArray_real_T *F
// Return Type  : void
//
static void activityFeatureExtractionAcc(const double sX[50], emxArray_real_T *F)
{
  emxArray_real_T *iPk;
  emxArray_real_T *idx;
  double dv0[49];
  double dv1[49];
  double dv2[49];
  int ixstart;
  emxArray_real_T *Min;
  emxArray_real_T *Imin;
  emxArray_real_T *b_idx;
  int i0;
  emxArray_real_T *Amplitude;
  emxArray_real_T *b_iPk;
  emxArray_real_T *Velocity;
  double y;
  double mtmp;
  int ix;
  boolean_T exitg2;
  double b_mtmp;
  boolean_T exitg1;
  emxArray_real_T *c_iPk;
  double T_average_peak_distance;
  double dv3[50];
  emxArray_real_T *d_iPk;
  int T_count_findpeaks;
  double T_findpeaks_distX;
  boolean_T x[50];
  double T_countmin_1;
  double T_countmin_2;
  double T_countmax;
  emxInit_real_T1(&iPk, 1);
  emxInit_real_T1(&idx, 1);

  // 0.175
  diff(sX, dv0);
  findpeaks(dv0, idx, iPk);
  diff(sX, dv1);
  for (ixstart = 0; ixstart < 49; ixstart++) {
    dv2[ixstart] = -dv1[ixstart];
  }

  emxInit_real_T1(&Min, 1);
  emxInit_real_T1(&Imin, 1);
  emxInit_real_T1(&b_idx, 1);
  findpeaks(dv2, Min, Imin);
  i0 = b_idx->size[0];
  b_idx->size[0] = idx->size[0];
  emxEnsureCapacity((emxArray__common *)b_idx, i0, (int)sizeof(double));
  ixstart = idx->size[0];
  for (i0 = 0; i0 < ixstart; i0++) {
    b_idx->data[i0] = idx->data[i0] - Min->data[i0];
  }

  emxInit_real_T1(&Amplitude, 1);
  emxInit_real_T1(&b_iPk, 1);
  b_abs(b_idx, Amplitude);
  i0 = b_iPk->size[0];
  b_iPk->size[0] = iPk->size[0];
  emxEnsureCapacity((emxArray__common *)b_iPk, i0, (int)sizeof(double));
  ixstart = iPk->size[0];
  emxFree_real_T(&b_idx);
  for (i0 = 0; i0 < ixstart; i0++) {
    b_iPk->data[i0] = iPk->data[i0] - Imin->data[i0];
  }

  emxInit_real_T(&Velocity, 2);
  mrdivide(Amplitude, b_iPk, Velocity);
  y = sX[0];
  emxFree_real_T(&b_iPk);
  for (ixstart = 0; ixstart < 49; ixstart++) {
    y += sX[ixstart + 1];
  }

  ixstart = 1;
  mtmp = sX[0];
  if (rtIsNaN(sX[0])) {
    ix = 2;
    exitg2 = false;
    while ((!exitg2) && (ix < 51)) {
      ixstart = ix;
      if (!rtIsNaN(sX[ix - 1])) {
        mtmp = sX[ix - 1];
        exitg2 = true;
      } else {
        ix++;
      }
    }
  }

  if (ixstart < 50) {
    while (ixstart + 1 < 51) {
      if (sX[ixstart] > mtmp) {
        mtmp = sX[ixstart];
      }

      ixstart++;
    }
  }

  ixstart = 1;
  b_mtmp = sX[0];
  if (rtIsNaN(sX[0])) {
    ix = 2;
    exitg1 = false;
    while ((!exitg1) && (ix < 51)) {
      ixstart = ix;
      if (!rtIsNaN(sX[ix - 1])) {
        b_mtmp = sX[ix - 1];
        exitg1 = true;
      } else {
        ix++;
      }
    }
  }

  if (ixstart < 50) {
    while (ixstart + 1 < 51) {
      if (sX[ixstart] < b_mtmp) {
        b_mtmp = sX[ixstart];
      }

      ixstart++;
    }
  }

  emxInit_real_T1(&c_iPk, 1);

  //  peaks = [];
  //  T_findpeaks_distX=[];
  T_average_peak_distance = 0.0;
  getAllPeaks(sX, iPk, Min, idx);
  removePeaksBelowMinPeakHeight(sX, iPk, 0.175);
  removePeaksBelowThreshold(sX, iPk, 0.0);
  combinePeaks(iPk, Min, c_iPk);
  c_findPeaksSeparatedByMoreThanM(c_iPk, idx);
  keepAtMostNpPeaks(idx);
  emxFree_real_T(&iPk);
  for (i0 = 0; i0 < 50; i0++) {
    dv3[i0] = 1.0 + (double)i0;
  }

  emxInit_real_T1(&d_iPk, 1);
  i0 = d_iPk->size[0];
  d_iPk->size[0] = idx->size[0];
  emxEnsureCapacity((emxArray__common *)d_iPk, i0, (int)sizeof(double));
  ixstart = idx->size[0];
  for (i0 = 0; i0 < ixstart; i0++) {
    d_iPk->data[i0] = c_iPk->data[(int)idx->data[i0] - 1];
  }

  emxFree_real_T(&c_iPk);
  assignOutputs(sX, dv3, d_iPk, Min, Imin);
  emxFree_real_T(&d_iPk);
  if (Min->size[0] == 0) {
    T_count_findpeaks = 0;
    T_findpeaks_distX = 0.0;
  } else {
    T_count_findpeaks = Min->size[0];
    i0 = idx->size[0];
    idx->size[0] = Min->size[0];
    emxEnsureCapacity((emxArray__common *)idx, i0, (int)sizeof(double));
    ixstart = Min->size[0];
    for (i0 = 0; i0 < ixstart; i0++) {
      idx->data[i0] = 0.0;
    }

    if (Min->size[0] > 1) {
      for (ixstart = 0; ixstart <= Min->size[0] - 2; ixstart++) {
        idx->data[ixstart] = Imin->data[1 + ixstart] - Imin->data[ixstart];
      }

      T_average_peak_distance = mean(idx);
      T_findpeaks_distX = Imin->data[Imin->size[0] - 1] - Imin->data[0];

      // TODO: TAKE AVG, NOT MAX-MIN
    } else {
      T_findpeaks_distX = 0.0;
    }
  }

  emxFree_real_T(&idx);
  emxFree_real_T(&Imin);
  emxFree_real_T(&Min);

  // %? Threshold Stuff:?
  for (ixstart = 0; ixstart < 50; ixstart++) {
    x[ixstart] = ((sX[ixstart] < 1.2) && (sX[ixstart] > 0.4));
  }

  T_countmin_1 = x[0];
  for (ixstart = 0; ixstart < 49; ixstart++) {
    T_countmin_1 += (double)x[ixstart + 1];
  }

  for (ixstart = 0; ixstart < 50; ixstart++) {
    x[ixstart] = (sX[ixstart] < 0.4);
  }

  T_countmin_2 = x[0];
  for (ixstart = 0; ixstart < 49; ixstart++) {
    T_countmin_2 += (double)x[ixstart + 1];
  }

  for (ixstart = 0; ixstart < 50; ixstart++) {
    x[ixstart] = (sX[ixstart] > 1.2);
  }

  T_countmax = x[0];
  for (ixstart = 0; ixstart < 49; ixstart++) {
    T_countmax += (double)x[ixstart + 1];
  }

  ixstart = !(Amplitude->size[0] == 0);
  if (!((Velocity->size[0] == 0) || (Velocity->size[1] == 0))) {
    ix = Velocity->size[1];
  } else {
    ix = 0;
  }

  i0 = F->size[0] * F->size[1];
  F->size[0] = 1;
  F->size[1] = (ixstart + ix) + 11;
  emxEnsureCapacity((emxArray__common *)F, i0, (int)sizeof(double));
  for (i0 = 0; i0 < ixstart; i0++) {
    F->data[F->size[0] * i0] = Amplitude->data[i0];
  }

  emxFree_real_T(&Amplitude);
  for (i0 = 0; i0 < ix; i0++) {
    F->data[F->size[0] * (i0 + ixstart)] = Velocity->data[i0];
  }

  emxFree_real_T(&Velocity);
  F->data[F->size[0] * (ixstart + ix)] = y;
  F->data[F->size[0] * ((ixstart + ix) + 1)] = b_std(sX);
  F->data[F->size[0] * ((ixstart + ix) + 2)] = mtmp;
  F->data[F->size[0] * ((ixstart + ix) + 3)] = b_mtmp;
  F->data[F->size[0] * ((ixstart + ix) + 4)] = trapz(sX);
  F->data[F->size[0] * ((ixstart + ix) + 5)] = T_count_findpeaks;
  F->data[F->size[0] * ((ixstart + ix) + 6)] = T_findpeaks_distX;
  F->data[F->size[0] * ((ixstart + ix) + 7)] = T_average_peak_distance;
  F->data[F->size[0] * ((ixstart + ix) + 8)] = T_countmax;
  F->data[F->size[0] * ((ixstart + ix) + 9)] = T_countmin_1;
  F->data[F->size[0] * ((ixstart + ix) + 10)] = T_countmin_2;

  //  hold off;
}

//
// featureExtraction Summary of this function goes here if I ever feel like
// writing one up.
//  sX = samples X input
// Arguments    : const double sX[50]
//                emxArray_real_T *F
// Return Type  : void
//
static void activityFeatureExtractionGyro(const double sX[50], emxArray_real_T
  *F)
{
  emxArray_real_T *iPk;
  emxArray_real_T *idx;
  double dv4[49];
  double dv5[49];
  double dv6[49];
  int ixstart;
  emxArray_real_T *Min;
  emxArray_real_T *Imin;
  emxArray_real_T *b_idx;
  int i6;
  emxArray_real_T *Amplitude;
  emxArray_real_T *b_iPk;
  emxArray_real_T *Velocity;
  double y;
  double mtmp;
  int ix;
  boolean_T exitg2;
  double b_mtmp;
  boolean_T exitg1;
  emxArray_real_T *c_iPk;
  double T_average_peak_distance;
  double dv7[50];
  emxArray_real_T *d_iPk;
  int T_count_findpeaks;
  double T_findpeaks_distX;
  boolean_T x[50];
  double T_countmin_1;
  double T_countmin_2;
  double T_countmax;
  emxInit_real_T1(&iPk, 1);
  emxInit_real_T1(&idx, 1);
  diff(sX, dv4);
  findpeaks(dv4, idx, iPk);
  diff(sX, dv5);
  for (ixstart = 0; ixstart < 49; ixstart++) {
    dv6[ixstart] = -dv5[ixstart];
  }

  emxInit_real_T1(&Min, 1);
  emxInit_real_T1(&Imin, 1);
  emxInit_real_T1(&b_idx, 1);
  findpeaks(dv6, Min, Imin);
  i6 = b_idx->size[0];
  b_idx->size[0] = idx->size[0];
  emxEnsureCapacity((emxArray__common *)b_idx, i6, (int)sizeof(double));
  ixstart = idx->size[0];
  for (i6 = 0; i6 < ixstart; i6++) {
    b_idx->data[i6] = idx->data[i6] - Min->data[i6];
  }

  emxInit_real_T1(&Amplitude, 1);
  emxInit_real_T1(&b_iPk, 1);
  b_abs(b_idx, Amplitude);
  i6 = b_iPk->size[0];
  b_iPk->size[0] = iPk->size[0];
  emxEnsureCapacity((emxArray__common *)b_iPk, i6, (int)sizeof(double));
  ixstart = iPk->size[0];
  emxFree_real_T(&b_idx);
  for (i6 = 0; i6 < ixstart; i6++) {
    b_iPk->data[i6] = iPk->data[i6] - Imin->data[i6];
  }

  emxInit_real_T(&Velocity, 2);
  mrdivide(Amplitude, b_iPk, Velocity);
  y = sX[0];
  emxFree_real_T(&b_iPk);
  for (ixstart = 0; ixstart < 49; ixstart++) {
    y += sX[ixstart + 1];
  }

  ixstart = 1;
  mtmp = sX[0];
  if (rtIsNaN(sX[0])) {
    ix = 2;
    exitg2 = false;
    while ((!exitg2) && (ix < 51)) {
      ixstart = ix;
      if (!rtIsNaN(sX[ix - 1])) {
        mtmp = sX[ix - 1];
        exitg2 = true;
      } else {
        ix++;
      }
    }
  }

  if (ixstart < 50) {
    while (ixstart + 1 < 51) {
      if (sX[ixstart] > mtmp) {
        mtmp = sX[ixstart];
      }

      ixstart++;
    }
  }

  ixstart = 1;
  b_mtmp = sX[0];
  if (rtIsNaN(sX[0])) {
    ix = 2;
    exitg1 = false;
    while ((!exitg1) && (ix < 51)) {
      ixstart = ix;
      if (!rtIsNaN(sX[ix - 1])) {
        b_mtmp = sX[ix - 1];
        exitg1 = true;
      } else {
        ix++;
      }
    }
  }

  if (ixstart < 50) {
    while (ixstart + 1 < 51) {
      if (sX[ixstart] < b_mtmp) {
        b_mtmp = sX[ixstart];
      }

      ixstart++;
    }
  }

  emxInit_real_T1(&c_iPk, 1);

  //  peaks = [];
  T_average_peak_distance = 0.0;

  //  T_findpeaks_distX=[];
  getAllPeaks(sX, iPk, Min, idx);
  removePeaksBelowMinPeakHeight(sX, iPk, 51.0);
  removePeaksBelowThreshold(sX, iPk, 0.0);
  combinePeaks(iPk, Min, c_iPk);
  c_findPeaksSeparatedByMoreThanM(c_iPk, idx);
  keepAtMostNpPeaks(idx);
  emxFree_real_T(&iPk);
  for (i6 = 0; i6 < 50; i6++) {
    dv7[i6] = 1.0 + (double)i6;
  }

  emxInit_real_T1(&d_iPk, 1);
  i6 = d_iPk->size[0];
  d_iPk->size[0] = idx->size[0];
  emxEnsureCapacity((emxArray__common *)d_iPk, i6, (int)sizeof(double));
  ixstart = idx->size[0];
  for (i6 = 0; i6 < ixstart; i6++) {
    d_iPk->data[i6] = c_iPk->data[(int)idx->data[i6] - 1];
  }

  emxFree_real_T(&c_iPk);
  assignOutputs(sX, dv7, d_iPk, Min, Imin);
  emxFree_real_T(&d_iPk);
  if (Min->size[0] == 0) {
    T_count_findpeaks = 0;
    T_findpeaks_distX = 0.0;
  } else {
    T_count_findpeaks = Min->size[0];
    i6 = idx->size[0];
    idx->size[0] = Min->size[0];
    emxEnsureCapacity((emxArray__common *)idx, i6, (int)sizeof(double));
    ixstart = Min->size[0];
    for (i6 = 0; i6 < ixstart; i6++) {
      idx->data[i6] = 0.0;
    }

    if (Min->size[0] > 1) {
      for (ixstart = 0; ixstart <= Min->size[0] - 2; ixstart++) {
        idx->data[ixstart] = Imin->data[1 + ixstart] - Imin->data[ixstart];
      }

      T_average_peak_distance = mean(idx);
      T_findpeaks_distX = Imin->data[Imin->size[0] - 1] - Imin->data[0];

      // TODO: TAKE AVG, NOT MAX-MIN
    } else {
      T_findpeaks_distX = 0.0;
    }
  }

  emxFree_real_T(&idx);
  emxFree_real_T(&Imin);
  emxFree_real_T(&Min);

  // %? Threshold Stuff:?
  for (ixstart = 0; ixstart < 50; ixstart++) {
    x[ixstart] = ((sX[ixstart] < 175.0) && (sX[ixstart] > 92.0));
  }

  T_countmin_1 = x[0];
  for (ixstart = 0; ixstart < 49; ixstart++) {
    T_countmin_1 += (double)x[ixstart + 1];
  }

  for (ixstart = 0; ixstart < 50; ixstart++) {
    x[ixstart] = (sX[ixstart] < 92.0);
  }

  T_countmin_2 = x[0];
  for (ixstart = 0; ixstart < 49; ixstart++) {
    T_countmin_2 += (double)x[ixstart + 1];
  }

  for (ixstart = 0; ixstart < 50; ixstart++) {
    x[ixstart] = (sX[ixstart] > 175.0);
  }

  T_countmax = x[0];
  for (ixstart = 0; ixstart < 49; ixstart++) {
    T_countmax += (double)x[ixstart + 1];
  }

  ixstart = !(Amplitude->size[0] == 0);
  if (!((Velocity->size[0] == 0) || (Velocity->size[1] == 0))) {
    ix = Velocity->size[1];
  } else {
    ix = 0;
  }

  i6 = F->size[0] * F->size[1];
  F->size[0] = 1;
  F->size[1] = (ixstart + ix) + 11;
  emxEnsureCapacity((emxArray__common *)F, i6, (int)sizeof(double));
  for (i6 = 0; i6 < ixstart; i6++) {
    F->data[F->size[0] * i6] = Amplitude->data[i6];
  }

  emxFree_real_T(&Amplitude);
  for (i6 = 0; i6 < ix; i6++) {
    F->data[F->size[0] * (i6 + ixstart)] = Velocity->data[i6];
  }

  emxFree_real_T(&Velocity);
  F->data[F->size[0] * (ixstart + ix)] = y;
  F->data[F->size[0] * ((ixstart + ix) + 1)] = b_std(sX);
  F->data[F->size[0] * ((ixstart + ix) + 2)] = mtmp;
  F->data[F->size[0] * ((ixstart + ix) + 3)] = b_mtmp;
  F->data[F->size[0] * ((ixstart + ix) + 4)] = trapz(sX);
  F->data[F->size[0] * ((ixstart + ix) + 5)] = T_count_findpeaks;
  F->data[F->size[0] * ((ixstart + ix) + 6)] = T_findpeaks_distX;
  F->data[F->size[0] * ((ixstart + ix) + 7)] = T_average_peak_distance;
  F->data[F->size[0] * ((ixstart + ix) + 8)] = T_countmax;
  F->data[F->size[0] * ((ixstart + ix) + 9)] = T_countmin_1;
  F->data[F->size[0] * ((ixstart + ix) + 10)] = T_countmin_2;

  //  hold off;
}

//
// Arguments    : const double y[50]
//                const double x[50]
//                const emxArray_real_T *iPk
//                emxArray_real_T *YpkOut
//                emxArray_real_T *XpkOut
// Return Type  : void
//
static void assignOutputs(const double y[50], const double x[50], const
  emxArray_real_T *iPk, emxArray_real_T *YpkOut, emxArray_real_T *XpkOut)
{
  int i5;
  int loop_ub;
  i5 = YpkOut->size[0];
  YpkOut->size[0] = iPk->size[0];
  emxEnsureCapacity((emxArray__common *)YpkOut, i5, (int)sizeof(double));
  loop_ub = iPk->size[0];
  for (i5 = 0; i5 < loop_ub; i5++) {
    YpkOut->data[i5] = y[(int)iPk->data[i5] - 1];
  }

  i5 = XpkOut->size[0];
  XpkOut->size[0] = iPk->size[0];
  emxEnsureCapacity((emxArray__common *)XpkOut, i5, (int)sizeof(double));
  loop_ub = iPk->size[0];
  for (i5 = 0; i5 < loop_ub; i5++) {
    XpkOut->data[i5] = x[(int)iPk->data[i5] - 1];
  }
}

//
// Arguments    : const emxArray_real_T *x
//                emxArray_real_T *y
// Return Type  : void
//
static void b_abs(const emxArray_real_T *x, emxArray_real_T *y)
{
  unsigned int unnamed_idx_0;
  int k;
  unnamed_idx_0 = (unsigned int)x->size[0];
  k = y->size[0];
  y->size[0] = (int)unnamed_idx_0;
  emxEnsureCapacity((emxArray__common *)y, k, (int)sizeof(double));
  for (k = 0; k + 1 <= x->size[0]; k++) {
    y->data[k] = std::abs(x->data[k]);
  }
}

//
// Arguments    : const emxArray_real_T *x
//                emxArray_real_T *y
// Return Type  : void
//
static void b_diff(const emxArray_real_T *x, emxArray_real_T *y)
{
  int iyLead;
  int orderForDim;
  emxArray_real_T *work;
  int ySize_idx_0;
  int m;
  double tmp1;
  int k;
  double tmp2;
  if (x->size[0] == 0) {
    iyLead = y->size[0];
    y->size[0] = 0;
    emxEnsureCapacity((emxArray__common *)y, iyLead, (int)sizeof(double));
  } else {
    if (x->size[0] - 1 <= 1) {
      orderForDim = x->size[0] - 1;
    } else {
      orderForDim = 1;
    }

    if (orderForDim < 1) {
      iyLead = y->size[0];
      y->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)y, iyLead, (int)sizeof(double));
    } else {
      emxInit_real_T1(&work, 1);
      ySize_idx_0 = x->size[0] - orderForDim;
      iyLead = work->size[0];
      work->size[0] = orderForDim;
      emxEnsureCapacity((emxArray__common *)work, iyLead, (int)sizeof(double));
      iyLead = y->size[0];
      y->size[0] = ySize_idx_0;
      emxEnsureCapacity((emxArray__common *)y, iyLead, (int)sizeof(double));
      if (!(y->size[0] == 0)) {
        ySize_idx_0 = 1;
        iyLead = 0;
        work->data[0] = x->data[0];
        if (orderForDim >= 2) {
          for (m = 1; m < orderForDim; m++) {
            tmp1 = x->data[ySize_idx_0];
            for (k = 0; k + 1 <= m; k++) {
              tmp2 = work->data[k];
              work->data[k] = tmp1;
              tmp1 -= tmp2;
            }

            work->data[m] = tmp1;
            ySize_idx_0++;
          }
        }

        for (m = orderForDim + 1; m <= x->size[0]; m++) {
          tmp1 = x->data[ySize_idx_0];
          for (k = 0; k + 1 <= orderForDim; k++) {
            tmp2 = work->data[k];
            work->data[k] = tmp1;
            tmp1 -= tmp2;
          }

          ySize_idx_0++;
          y->data[iyLead] = tmp1;
          iyLead++;
        }
      }

      emxFree_real_T(&work);
    }
  }
}

//
// Arguments    : const double yTemp[50]
//                emxArray_real_T *iPk
//                emxArray_real_T *iInflect
// Return Type  : void
//
static void b_findLocalMaxima(const double yTemp[50], emxArray_real_T *iPk,
  emxArray_real_T *iInflect)
{
  double b_yTemp[52];
  boolean_T yFinite[52];
  int ii;
  boolean_T x[51];
  emxArray_int32_T *b_ii;
  int idx;
  int i4;
  boolean_T exitg3;
  emxArray_int32_T *r2;
  boolean_T guard3 = false;
  emxArray_real_T *iTemp;
  emxArray_real_T *c_yTemp;
  emxArray_real_T *s;
  emxArray_boolean_T *b_x;
  emxArray_real_T *r3;
  int nx;
  boolean_T exitg2;
  boolean_T guard2 = false;
  emxArray_int32_T *c_ii;
  boolean_T exitg1;
  boolean_T guard1 = false;
  b_yTemp[0] = rtNaN;
  memcpy(&b_yTemp[1], &yTemp[0], 50U * sizeof(double));
  b_yTemp[51] = rtNaN;
  for (ii = 0; ii < 52; ii++) {
    yFinite[ii] = !rtIsNaN(b_yTemp[ii]);
  }

  for (ii = 0; ii < 51; ii++) {
    x[ii] = ((b_yTemp[ii] != b_yTemp[ii + 1]) && (yFinite[ii] || yFinite[ii + 1]));
  }

  emxInit_int32_T(&b_ii, 1);
  idx = 0;
  i4 = b_ii->size[0];
  b_ii->size[0] = 51;
  emxEnsureCapacity((emxArray__common *)b_ii, i4, (int)sizeof(int));
  ii = 1;
  exitg3 = false;
  while ((!exitg3) && (ii < 52)) {
    guard3 = false;
    if (x[ii - 1]) {
      idx++;
      b_ii->data[idx - 1] = ii;
      if (idx >= 51) {
        exitg3 = true;
      } else {
        guard3 = true;
      }
    } else {
      guard3 = true;
    }

    if (guard3) {
      ii++;
    }
  }

  emxInit_int32_T(&r2, 1);
  i4 = b_ii->size[0];
  if (1 > idx) {
    b_ii->size[0] = 0;
  } else {
    b_ii->size[0] = idx;
  }

  emxEnsureCapacity((emxArray__common *)b_ii, i4, (int)sizeof(int));
  i4 = r2->size[0];
  r2->size[0] = 1 + b_ii->size[0];
  emxEnsureCapacity((emxArray__common *)r2, i4, (int)sizeof(int));
  r2->data[0] = 1;
  ii = b_ii->size[0];
  for (i4 = 0; i4 < ii; i4++) {
    r2->data[i4 + 1] = b_ii->data[i4] + 1;
  }

  emxInit_real_T1(&iTemp, 1);
  i4 = iTemp->size[0];
  iTemp->size[0] = r2->size[0];
  emxEnsureCapacity((emxArray__common *)iTemp, i4, (int)sizeof(double));
  ii = r2->size[0];
  for (i4 = 0; i4 < ii; i4++) {
    iTemp->data[i4] = 1.0 + (double)(r2->data[i4] - 1);
  }

  emxFree_int32_T(&r2);
  emxInit_real_T1(&c_yTemp, 1);
  i4 = c_yTemp->size[0];
  c_yTemp->size[0] = iTemp->size[0];
  emxEnsureCapacity((emxArray__common *)c_yTemp, i4, (int)sizeof(double));
  ii = iTemp->size[0];
  for (i4 = 0; i4 < ii; i4++) {
    c_yTemp->data[i4] = b_yTemp[(int)iTemp->data[i4] - 1];
  }

  emxInit_real_T1(&s, 1);
  emxInit_boolean_T(&b_x, 1);
  emxInit_real_T1(&r3, 1);
  b_diff(c_yTemp, s);
  b_sign(s);
  b_diff(s, r3);
  i4 = b_x->size[0];
  b_x->size[0] = r3->size[0];
  emxEnsureCapacity((emxArray__common *)b_x, i4, (int)sizeof(boolean_T));
  ii = r3->size[0];
  emxFree_real_T(&c_yTemp);
  for (i4 = 0; i4 < ii; i4++) {
    b_x->data[i4] = (r3->data[i4] < 0.0);
  }

  emxFree_real_T(&r3);
  nx = b_x->size[0];
  idx = 0;
  i4 = b_ii->size[0];
  b_ii->size[0] = b_x->size[0];
  emxEnsureCapacity((emxArray__common *)b_ii, i4, (int)sizeof(int));
  ii = 1;
  exitg2 = false;
  while ((!exitg2) && (ii <= nx)) {
    guard2 = false;
    if (b_x->data[ii - 1]) {
      idx++;
      b_ii->data[idx - 1] = ii;
      if (idx >= nx) {
        exitg2 = true;
      } else {
        guard2 = true;
      }
    } else {
      guard2 = true;
    }

    if (guard2) {
      ii++;
    }
  }

  if (b_x->size[0] == 1) {
    if (idx == 0) {
      i4 = b_ii->size[0];
      b_ii->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)b_ii, i4, (int)sizeof(int));
    }
  } else {
    i4 = b_ii->size[0];
    if (1 > idx) {
      b_ii->size[0] = 0;
    } else {
      b_ii->size[0] = idx;
    }

    emxEnsureCapacity((emxArray__common *)b_ii, i4, (int)sizeof(int));
  }

  if (1.0 > (double)s->size[0] - 1.0) {
    ii = 0;
  } else {
    ii = (int)((double)s->size[0] - 1.0);
  }

  if (2 > s->size[0]) {
    i4 = 0;
  } else {
    i4 = 1;
  }

  idx = b_x->size[0];
  b_x->size[0] = ii;
  emxEnsureCapacity((emxArray__common *)b_x, idx, (int)sizeof(boolean_T));
  for (idx = 0; idx < ii; idx++) {
    b_x->data[idx] = (s->data[idx] != s->data[i4 + idx]);
  }

  emxFree_real_T(&s);
  emxInit_int32_T(&c_ii, 1);
  nx = b_x->size[0];
  idx = 0;
  i4 = c_ii->size[0];
  c_ii->size[0] = b_x->size[0];
  emxEnsureCapacity((emxArray__common *)c_ii, i4, (int)sizeof(int));
  ii = 1;
  exitg1 = false;
  while ((!exitg1) && (ii <= nx)) {
    guard1 = false;
    if (b_x->data[ii - 1]) {
      idx++;
      c_ii->data[idx - 1] = ii;
      if (idx >= nx) {
        exitg1 = true;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1) {
      ii++;
    }
  }

  if (b_x->size[0] == 1) {
    if (idx == 0) {
      i4 = c_ii->size[0];
      c_ii->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)c_ii, i4, (int)sizeof(int));
    }
  } else {
    i4 = c_ii->size[0];
    if (1 > idx) {
      c_ii->size[0] = 0;
    } else {
      c_ii->size[0] = idx;
    }

    emxEnsureCapacity((emxArray__common *)c_ii, i4, (int)sizeof(int));
  }

  emxFree_boolean_T(&b_x);
  i4 = iInflect->size[0];
  iInflect->size[0] = c_ii->size[0];
  emxEnsureCapacity((emxArray__common *)iInflect, i4, (int)sizeof(double));
  ii = c_ii->size[0];
  for (i4 = 0; i4 < ii; i4++) {
    iInflect->data[i4] = iTemp->data[c_ii->data[i4]] - 1.0;
  }

  emxFree_int32_T(&c_ii);
  i4 = iPk->size[0];
  iPk->size[0] = b_ii->size[0];
  emxEnsureCapacity((emxArray__common *)iPk, i4, (int)sizeof(double));
  ii = b_ii->size[0];
  for (i4 = 0; i4 < ii; i4++) {
    iPk->data[i4] = iTemp->data[b_ii->data[i4]] - 1.0;
  }

  emxFree_int32_T(&b_ii);
  emxFree_real_T(&iTemp);
}

//
// Arguments    : int idx[245]
//                double x[245]
//                int offset
//                int np
//                int nq
//                int iwork[245]
//                double xwork[245]
// Return Type  : void
//
static void b_merge(int idx[245], double x[245], int offset, int np, int nq, int
                    iwork[245], double xwork[245])
{
  int n;
  int qend;
  int p;
  int iout;
  int exitg1;
  if ((np == 0) || (nq == 0)) {
  } else {
    n = np + nq;
    for (qend = 0; qend + 1 <= n; qend++) {
      iwork[qend] = idx[offset + qend];
      xwork[qend] = x[offset + qend];
    }

    p = 0;
    n = np;
    qend = np + nq;
    iout = offset - 1;
    do {
      exitg1 = 0;
      iout++;
      if (xwork[p] <= xwork[n]) {
        idx[iout] = iwork[p];
        x[iout] = xwork[p];
        if (p + 1 < np) {
          p++;
        } else {
          exitg1 = 1;
        }
      } else {
        idx[iout] = iwork[n];
        x[iout] = xwork[n];
        if (n + 1 < qend) {
          n++;
        } else {
          n = iout - p;
          while (p + 1 <= np) {
            idx[(n + p) + 1] = iwork[p];
            x[(n + p) + 1] = xwork[p];
            p++;
          }

          exitg1 = 1;
        }
      }
    } while (exitg1 == 0);
  }
}

//
// Arguments    : int idx[245]
//                double x[245]
//                int offset
//                int n
//                int preSortLevel
//                int iwork[245]
//                double xwork[245]
// Return Type  : void
//
static void b_merge_block(int idx[245], double x[245], int offset, int n, int
  preSortLevel, int iwork[245], double xwork[245])
{
  int nPairs;
  int bLen;
  int tailOffset;
  int nTail;
  nPairs = n >> preSortLevel;
  bLen = 1 << preSortLevel;
  while (nPairs > 1) {
    if ((nPairs & 1) != 0) {
      nPairs--;
      tailOffset = bLen * nPairs;
      nTail = n - tailOffset;
      if (nTail > bLen) {
        b_merge(idx, x, offset + tailOffset, bLen, nTail - bLen, iwork, xwork);
      }
    }

    tailOffset = bLen << 1;
    nPairs >>= 1;
    for (nTail = 1; nTail <= nPairs; nTail++) {
      b_merge(idx, x, offset + (nTail - 1) * tailOffset, bLen, bLen, iwork,
              xwork);
    }

    bLen = tailOffset;
  }

  if (n > bLen) {
    b_merge(idx, x, offset, bLen, n - bLen, iwork, xwork);
  }
}

//
// Arguments    : emxArray_real_T *x
// Return Type  : void
//
static void b_sign(emxArray_real_T *x)
{
  int nx;
  int k;
  double b_x;
  nx = x->size[0];
  for (k = 0; k + 1 <= nx; k++) {
    if (x->data[k] < 0.0) {
      b_x = -1.0;
    } else if (x->data[k] > 0.0) {
      b_x = 1.0;
    } else if (x->data[k] == 0.0) {
      b_x = 0.0;
    } else {
      b_x = x->data[k];
    }

    x->data[k] = b_x;
  }
}

//
// Arguments    : emxArray_real_T *x
//                int dim
//                emxArray_int32_T *idx
// Return Type  : void
//
static void b_sort(emxArray_real_T *x, int dim, emxArray_int32_T *idx)
{
  int i8;
  emxArray_real_T *vwork;
  int k;
  unsigned int unnamed_idx_0;
  int vstride;
  emxArray_int32_T *iidx;
  int j;
  if (dim <= 1) {
    i8 = x->size[0];
  } else {
    i8 = 1;
  }

  emxInit_real_T1(&vwork, 1);
  k = vwork->size[0];
  vwork->size[0] = i8;
  emxEnsureCapacity((emxArray__common *)vwork, k, (int)sizeof(double));
  unnamed_idx_0 = (unsigned int)x->size[0];
  k = idx->size[0];
  idx->size[0] = (int)unnamed_idx_0;
  emxEnsureCapacity((emxArray__common *)idx, k, (int)sizeof(int));
  if (dim > 2) {
    vstride = x->size[0];
  } else {
    vstride = 1;
    k = 1;
    while (k <= dim - 1) {
      k = x->size[0];
      vstride *= k;
      k = 2;
    }
  }

  emxInit_int32_T(&iidx, 1);
  for (j = 0; j + 1 <= vstride; j++) {
    for (k = 0; k + 1 <= i8; k++) {
      vwork->data[k] = x->data[j + k * vstride];
    }

    sortIdx(vwork, iidx);
    for (k = 0; k + 1 <= i8; k++) {
      x->data[j + k * vstride] = vwork->data[k];
      idx->data[j + k * vstride] = iidx->data[k];
    }
  }

  emxFree_int32_T(&iidx);
  emxFree_real_T(&vwork);
}

//
// Arguments    : const double varargin_1[50]
// Return Type  : double
//
static double b_std(const double varargin_1[50])
{
  double y;
  int ix;
  double xbar;
  int k;
  double r;
  ix = 0;
  xbar = varargin_1[0];
  for (k = 0; k < 49; k++) {
    ix++;
    xbar += varargin_1[ix];
  }

  xbar /= 50.0;
  ix = 0;
  r = varargin_1[0] - xbar;
  y = r * r;
  for (k = 0; k < 49; k++) {
    ix++;
    r = varargin_1[ix] - xbar;
    y += r * r;
  }

  y /= 49.0;
  return std::sqrt(y);
}

//
// Arguments    : const emxArray_real_T *iPk
//                emxArray_real_T *idx
// Return Type  : void
//
static void c_findPeaksSeparatedByMoreThanM(const emxArray_real_T *iPk,
  emxArray_real_T *idx)
{
  int ndbl;
  int apnd;
  int cdiff;
  int absb;
  emxArray_real_T *y;
  if (iPk->size[0] < 1) {
    ndbl = 0;
    apnd = 0;
  } else {
    ndbl = (int)std::floor(((double)iPk->size[0] - 1.0) + 0.5);
    apnd = ndbl + 1;
    cdiff = (ndbl - iPk->size[0]) + 1;
    absb = iPk->size[0];
    if (std::abs((double)cdiff) < 4.4408920985006262E-16 * (double)absb) {
      ndbl++;
      apnd = iPk->size[0];
    } else if (cdiff > 0) {
      apnd = ndbl;
    } else {
      ndbl++;
    }
  }

  emxInit_real_T(&y, 2);
  absb = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = ndbl;
  emxEnsureCapacity((emxArray__common *)y, absb, (int)sizeof(double));
  if (ndbl > 0) {
    y->data[0] = 1.0;
    if (ndbl > 1) {
      y->data[ndbl - 1] = apnd;
      cdiff = (ndbl - 1) / 2;
      for (absb = 1; absb < cdiff; absb++) {
        y->data[absb] = 1.0 + (double)absb;
        y->data[(ndbl - absb) - 1] = apnd - absb;
      }

      if (cdiff << 1 == ndbl - 1) {
        y->data[cdiff] = (1.0 + (double)apnd) / 2.0;
      } else {
        y->data[cdiff] = 1.0 + (double)cdiff;
        y->data[cdiff + 1] = apnd - cdiff;
      }
    }
  }

  absb = idx->size[0];
  idx->size[0] = y->size[1];
  emxEnsureCapacity((emxArray__common *)idx, absb, (int)sizeof(double));
  cdiff = y->size[1];
  for (absb = 0; absb < cdiff; absb++) {
    idx->data[absb] = y->data[y->size[0] * absb];
  }

  emxFree_real_T(&y);
}

//
// Arguments    : double x[245]
//                int idx[245]
// Return Type  : void
//
static void c_sort(double x[245], int idx[245])
{
  double x4[4];
  unsigned char idx4[4];
  int m;
  double xwork[245];
  int nNaNs;
  int ib;
  int k;
  signed char perm[4];
  int i2;
  int iwork[245];
  int i3;
  int i4;
  memset(&idx[0], 0, 245U * sizeof(int));
  for (m = 0; m < 4; m++) {
    x4[m] = 0.0;
    idx4[m] = 0;
  }

  memset(&xwork[0], 0, 245U * sizeof(double));
  nNaNs = -244;
  ib = 0;
  for (k = 0; k < 245; k++) {
    if (rtIsNaN(x[k])) {
      idx[-nNaNs] = k + 1;
      xwork[-nNaNs] = x[k];
      nNaNs++;
    } else {
      ib++;
      idx4[ib - 1] = (unsigned char)(k + 1);
      x4[ib - 1] = x[k];
      if (ib == 4) {
        ib = (k - nNaNs) - 247;
        if (x4[0] <= x4[1]) {
          m = 1;
          i2 = 2;
        } else {
          m = 2;
          i2 = 1;
        }

        if (x4[2] <= x4[3]) {
          i3 = 3;
          i4 = 4;
        } else {
          i3 = 4;
          i4 = 3;
        }

        if (x4[m - 1] <= x4[i3 - 1]) {
          if (x4[i2 - 1] <= x4[i3 - 1]) {
            perm[0] = (signed char)m;
            perm[1] = (signed char)i2;
            perm[2] = (signed char)i3;
            perm[3] = (signed char)i4;
          } else if (x4[i2 - 1] <= x4[i4 - 1]) {
            perm[0] = (signed char)m;
            perm[1] = (signed char)i3;
            perm[2] = (signed char)i2;
            perm[3] = (signed char)i4;
          } else {
            perm[0] = (signed char)m;
            perm[1] = (signed char)i3;
            perm[2] = (signed char)i4;
            perm[3] = (signed char)i2;
          }
        } else if (x4[m - 1] <= x4[i4 - 1]) {
          if (x4[i2 - 1] <= x4[i4 - 1]) {
            perm[0] = (signed char)i3;
            perm[1] = (signed char)m;
            perm[2] = (signed char)i2;
            perm[3] = (signed char)i4;
          } else {
            perm[0] = (signed char)i3;
            perm[1] = (signed char)m;
            perm[2] = (signed char)i4;
            perm[3] = (signed char)i2;
          }
        } else {
          perm[0] = (signed char)i3;
          perm[1] = (signed char)i4;
          perm[2] = (signed char)m;
          perm[3] = (signed char)i2;
        }

        idx[ib] = idx4[perm[0] - 1];
        idx[ib + 1] = idx4[perm[1] - 1];
        idx[ib + 2] = idx4[perm[2] - 1];
        idx[ib + 3] = idx4[perm[3] - 1];
        x[ib] = x4[perm[0] - 1];
        x[ib + 1] = x4[perm[1] - 1];
        x[ib + 2] = x4[perm[2] - 1];
        x[ib + 3] = x4[perm[3] - 1];
        ib = 0;
      }
    }
  }

  if (ib > 0) {
    for (m = 0; m < 4; m++) {
      perm[m] = 0;
    }

    if (ib == 1) {
      perm[0] = 1;
    } else if (ib == 2) {
      if (x4[0] <= x4[1]) {
        perm[0] = 1;
        perm[1] = 2;
      } else {
        perm[0] = 2;
        perm[1] = 1;
      }
    } else if (x4[0] <= x4[1]) {
      if (x4[1] <= x4[2]) {
        perm[0] = 1;
        perm[1] = 2;
        perm[2] = 3;
      } else if (x4[0] <= x4[2]) {
        perm[0] = 1;
        perm[1] = 3;
        perm[2] = 2;
      } else {
        perm[0] = 3;
        perm[1] = 1;
        perm[2] = 2;
      }
    } else if (x4[0] <= x4[2]) {
      perm[0] = 2;
      perm[1] = 1;
      perm[2] = 3;
    } else if (x4[1] <= x4[2]) {
      perm[0] = 2;
      perm[1] = 3;
      perm[2] = 1;
    } else {
      perm[0] = 3;
      perm[1] = 2;
      perm[2] = 1;
    }

    for (k = 1; k <= ib; k++) {
      idx[(k - nNaNs) - ib] = idx4[perm[k - 1] - 1];
      x[(k - nNaNs) - ib] = x4[perm[k - 1] - 1];
    }
  }

  m = (nNaNs + 244) >> 1;
  for (k = 1; k <= m; k++) {
    ib = idx[k - nNaNs];
    idx[k - nNaNs] = idx[245 - k];
    idx[245 - k] = ib;
    x[k - nNaNs] = xwork[245 - k];
    x[245 - k] = xwork[k - nNaNs];
  }

  if (((nNaNs + 244) & 1) != 0) {
    x[(m - nNaNs) + 1] = xwork[(m - nNaNs) + 1];
  }

  if (1 - nNaNs > 1) {
    memset(&iwork[0], 0, 245U * sizeof(int));
    b_merge_block(idx, x, 0, 1 - nNaNs, 2, iwork, xwork);
  }
}

//
// Arguments    : const emxArray_real_T *iPk
//                const emxArray_real_T *iInf
//                emxArray_real_T *iPkOut
// Return Type  : void
//
static void combinePeaks(const emxArray_real_T *iPk, const emxArray_real_T *iInf,
  emxArray_real_T *iPkOut)
{
  emxArray_int32_T *ia;
  emxArray_int32_T *ib;
  emxInit_int32_T(&ia, 1);
  emxInit_int32_T(&ib, 1);
  do_vectors(iPk, iInf, iPkOut, ia, ib);
  emxFree_int32_T(&ib);
  emxFree_int32_T(&ia);
}

//
// Arguments    : const double x[50]
//                double y[49]
// Return Type  : void
//
static void diff(const double x[50], double y[49])
{
  int ixLead;
  int iyLead;
  double work;
  int m;
  double tmp2;
  ixLead = 1;
  iyLead = 0;
  work = x[0];
  for (m = 0; m < 49; m++) {
    tmp2 = work;
    work = x[ixLead];
    tmp2 = x[ixLead] - tmp2;
    ixLead++;
    y[iyLead] = tmp2;
    iyLead++;
  }
}

//
// Arguments    : const emxArray_real_T *a
//                const emxArray_real_T *b
//                emxArray_real_T *c
//                emxArray_int32_T *ia
//                emxArray_int32_T *ib
// Return Type  : void
//
static void do_vectors(const emxArray_real_T *a, const emxArray_real_T *b,
  emxArray_real_T *c, emxArray_int32_T *ia, emxArray_int32_T *ib)
{
  int na;
  int nb;
  int ncmax;
  int ibfirst;
  int nc;
  int nia;
  int nib;
  int iafirst;
  int ialast;
  int iblast;
  int b_ialast;
  double ak;
  int b_iblast;
  double bk;
  boolean_T p;
  emxArray_int32_T *b_ia;
  emxArray_int32_T *b_ib;
  emxArray_real_T *b_c;
  na = a->size[0];
  nb = b->size[0];
  ncmax = a->size[0] + b->size[0];
  ibfirst = c->size[0];
  c->size[0] = ncmax;
  emxEnsureCapacity((emxArray__common *)c, ibfirst, (int)sizeof(double));
  ibfirst = ia->size[0];
  ia->size[0] = a->size[0];
  emxEnsureCapacity((emxArray__common *)ia, ibfirst, (int)sizeof(int));
  ibfirst = ib->size[0];
  ib->size[0] = b->size[0];
  emxEnsureCapacity((emxArray__common *)ib, ibfirst, (int)sizeof(int));
  nc = -1;
  nia = -1;
  nib = 0;
  iafirst = 1;
  ialast = 1;
  ibfirst = 0;
  iblast = 1;
  while ((ialast <= na) && (iblast <= nb)) {
    b_ialast = ialast;
    ak = skip_to_last_equal_value(&b_ialast, a);
    ialast = b_ialast;
    b_iblast = iblast;
    bk = skip_to_last_equal_value(&b_iblast, b);
    iblast = b_iblast;
    if ((std::abs(bk - ak) < eps(bk / 2.0)) || (rtIsInf(ak) && rtIsInf(bk) &&
         ((ak > 0.0) == (bk > 0.0)))) {
      p = true;
    } else {
      p = false;
    }

    if (p) {
      nc++;
      c->data[nc] = ak;
      nia++;
      ia->data[nia] = iafirst;
      ialast = b_ialast + 1;
      iafirst = b_ialast + 1;
      iblast = b_iblast + 1;
      ibfirst = b_iblast;
    } else {
      if ((ak < bk) || rtIsNaN(bk)) {
        p = true;
      } else {
        p = false;
      }

      if (p) {
        nc++;
        nia++;
        c->data[nc] = ak;
        ia->data[nia] = iafirst;
        ialast = b_ialast + 1;
        iafirst = b_ialast + 1;
      } else {
        nc++;
        nib++;
        c->data[nc] = bk;
        ib->data[nib - 1] = ibfirst + 1;
        iblast = b_iblast + 1;
        ibfirst = b_iblast;
      }
    }
  }

  while (ialast <= na) {
    iafirst = ialast;
    ak = skip_to_last_equal_value(&iafirst, a);
    nc++;
    nia++;
    c->data[nc] = ak;
    ia->data[nia] = ialast;
    ialast = iafirst + 1;
  }

  while (iblast <= nb) {
    iafirst = iblast;
    bk = skip_to_last_equal_value(&iafirst, b);
    nc++;
    nib++;
    c->data[nc] = bk;
    ib->data[nib - 1] = iblast;
    iblast = iafirst + 1;
  }

  if (a->size[0] > 0) {
    if (1 > nia + 1) {
      iafirst = -1;
    } else {
      iafirst = nia;
    }

    emxInit_int32_T(&b_ia, 1);
    ibfirst = b_ia->size[0];
    b_ia->size[0] = iafirst + 1;
    emxEnsureCapacity((emxArray__common *)b_ia, ibfirst, (int)sizeof(int));
    for (ibfirst = 0; ibfirst <= iafirst; ibfirst++) {
      b_ia->data[ibfirst] = ia->data[ibfirst];
    }

    ibfirst = ia->size[0];
    ia->size[0] = b_ia->size[0];
    emxEnsureCapacity((emxArray__common *)ia, ibfirst, (int)sizeof(int));
    iafirst = b_ia->size[0];
    for (ibfirst = 0; ibfirst < iafirst; ibfirst++) {
      ia->data[ibfirst] = b_ia->data[ibfirst];
    }

    emxFree_int32_T(&b_ia);
  }

  if (b->size[0] > 0) {
    if (1 > nib) {
      iafirst = 0;
    } else {
      iafirst = nib;
    }

    emxInit_int32_T(&b_ib, 1);
    ibfirst = b_ib->size[0];
    b_ib->size[0] = iafirst;
    emxEnsureCapacity((emxArray__common *)b_ib, ibfirst, (int)sizeof(int));
    for (ibfirst = 0; ibfirst < iafirst; ibfirst++) {
      b_ib->data[ibfirst] = ib->data[ibfirst];
    }

    ibfirst = ib->size[0];
    ib->size[0] = b_ib->size[0];
    emxEnsureCapacity((emxArray__common *)ib, ibfirst, (int)sizeof(int));
    iafirst = b_ib->size[0];
    for (ibfirst = 0; ibfirst < iafirst; ibfirst++) {
      ib->data[ibfirst] = b_ib->data[ibfirst];
    }

    emxFree_int32_T(&b_ib);
  }

  if (ncmax > 0) {
    if (1 > nc + 1) {
      iafirst = -1;
    } else {
      iafirst = nc;
    }

    emxInit_real_T1(&b_c, 1);
    ibfirst = b_c->size[0];
    b_c->size[0] = iafirst + 1;
    emxEnsureCapacity((emxArray__common *)b_c, ibfirst, (int)sizeof(double));
    for (ibfirst = 0; ibfirst <= iafirst; ibfirst++) {
      b_c->data[ibfirst] = c->data[ibfirst];
    }

    ibfirst = c->size[0];
    c->size[0] = b_c->size[0];
    emxEnsureCapacity((emxArray__common *)c, ibfirst, (int)sizeof(double));
    iafirst = b_c->size[0];
    for (ibfirst = 0; ibfirst < iafirst; ibfirst++) {
      c->data[ibfirst] = b_c->data[ibfirst];
    }

    emxFree_real_T(&b_c);
  }
}

//
// Arguments    : emxArray__common *emxArray
//                int oldNumel
//                int elementSize
// Return Type  : void
//
static void emxEnsureCapacity(emxArray__common *emxArray, int oldNumel, int
  elementSize)
{
  int newNumel;
  int i;
  void *newData;
  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }

  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }

    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i <<= 1;
      }
    }

    newData = calloc((unsigned int)i, (unsigned int)elementSize);
    if (emxArray->data != NULL) {
      memcpy(newData, emxArray->data, (unsigned int)(elementSize * oldNumel));
      if (emxArray->canFreeData) {
        free(emxArray->data);
      }
    }

    emxArray->data = newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

//
// Arguments    : emxArray_boolean_T **pEmxArray
// Return Type  : void
//
static void emxFree_boolean_T(emxArray_boolean_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_boolean_T *)NULL) {
    if (((*pEmxArray)->data != (boolean_T *)NULL) && (*pEmxArray)->canFreeData)
    {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_boolean_T *)NULL;
  }
}

//
// Arguments    : emxArray_int32_T **pEmxArray
// Return Type  : void
//
static void emxFree_int32_T(emxArray_int32_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_int32_T *)NULL) {
    if (((*pEmxArray)->data != (int *)NULL) && (*pEmxArray)->canFreeData) {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_int32_T *)NULL;
  }
}

//
// Arguments    : emxArray_real_T **pEmxArray
// Return Type  : void
//
static void emxFree_real_T(emxArray_real_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_real_T *)NULL) {
    if (((*pEmxArray)->data != (double *)NULL) && (*pEmxArray)->canFreeData) {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_real_T *)NULL;
  }
}

//
// Arguments    : emxArray_boolean_T **pEmxArray
//                int numDimensions
// Return Type  : void
//
static void emxInit_boolean_T(emxArray_boolean_T **pEmxArray, int numDimensions)
{
  emxArray_boolean_T *emxArray;
  int i;
  *pEmxArray = (emxArray_boolean_T *)malloc(sizeof(emxArray_boolean_T));
  emxArray = *pEmxArray;
  emxArray->data = (boolean_T *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

//
// Arguments    : emxArray_int32_T **pEmxArray
//                int numDimensions
// Return Type  : void
//
static void emxInit_int32_T(emxArray_int32_T **pEmxArray, int numDimensions)
{
  emxArray_int32_T *emxArray;
  int i;
  *pEmxArray = (emxArray_int32_T *)malloc(sizeof(emxArray_int32_T));
  emxArray = *pEmxArray;
  emxArray->data = (int *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

//
// Arguments    : emxArray_int32_T **pEmxArray
//                int numDimensions
// Return Type  : void
//
static void emxInit_int32_T1(emxArray_int32_T **pEmxArray, int numDimensions)
{
  emxArray_int32_T *emxArray;
  int i;
  *pEmxArray = (emxArray_int32_T *)malloc(sizeof(emxArray_int32_T));
  emxArray = *pEmxArray;
  emxArray->data = (int *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

//
// Arguments    : emxArray_real_T **pEmxArray
//                int numDimensions
// Return Type  : void
//
static void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions)
{
  emxArray_real_T *emxArray;
  int i;
  *pEmxArray = (emxArray_real_T *)malloc(sizeof(emxArray_real_T));
  emxArray = *pEmxArray;
  emxArray->data = (double *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

//
// Arguments    : emxArray_real_T **pEmxArray
//                int numDimensions
// Return Type  : void
//
static void emxInit_real_T1(emxArray_real_T **pEmxArray, int numDimensions)
{
  emxArray_real_T *emxArray;
  int i;
  *pEmxArray = (emxArray_real_T *)malloc(sizeof(emxArray_real_T));
  emxArray = *pEmxArray;
  emxArray->data = (double *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numDimensions));
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

//
// Arguments    : double x
// Return Type  : double
//
static double eps(double x)
{
  double r;
  double absxk;
  int exponent;
  absxk = std::abs(x);
  if ((!rtIsInf(absxk)) && (!rtIsNaN(absxk))) {
    if (absxk <= 2.2250738585072014E-308) {
      r = 4.94065645841247E-324;
    } else {
      frexp(absxk, &exponent);
      r = std::ldexp(1.0, exponent - 53);
    }
  } else {
    r = rtNaN;
  }

  return r;
}

//
// Arguments    : const double yTemp[49]
//                emxArray_real_T *iPk
//                emxArray_real_T *iInflect
// Return Type  : void
//
static void findLocalMaxima(const double yTemp[49], emxArray_real_T *iPk,
  emxArray_real_T *iInflect)
{
  double b_yTemp[51];
  boolean_T yFinite[51];
  int ii;
  boolean_T x[50];
  emxArray_int32_T *b_ii;
  int idx;
  int i1;
  boolean_T exitg3;
  emxArray_int32_T *r0;
  boolean_T guard3 = false;
  emxArray_real_T *iTemp;
  emxArray_real_T *c_yTemp;
  emxArray_real_T *s;
  emxArray_boolean_T *b_x;
  emxArray_real_T *r1;
  int nx;
  boolean_T exitg2;
  boolean_T guard2 = false;
  emxArray_int32_T *c_ii;
  boolean_T exitg1;
  boolean_T guard1 = false;
  b_yTemp[0] = rtNaN;
  memcpy(&b_yTemp[1], &yTemp[0], 49U * sizeof(double));
  b_yTemp[50] = rtNaN;
  for (ii = 0; ii < 51; ii++) {
    yFinite[ii] = !rtIsNaN(b_yTemp[ii]);
  }

  for (ii = 0; ii < 50; ii++) {
    x[ii] = ((b_yTemp[ii] != b_yTemp[ii + 1]) && (yFinite[ii] || yFinite[ii + 1]));
  }

  emxInit_int32_T(&b_ii, 1);
  idx = 0;
  i1 = b_ii->size[0];
  b_ii->size[0] = 50;
  emxEnsureCapacity((emxArray__common *)b_ii, i1, (int)sizeof(int));
  ii = 1;
  exitg3 = false;
  while ((!exitg3) && (ii < 51)) {
    guard3 = false;
    if (x[ii - 1]) {
      idx++;
      b_ii->data[idx - 1] = ii;
      if (idx >= 50) {
        exitg3 = true;
      } else {
        guard3 = true;
      }
    } else {
      guard3 = true;
    }

    if (guard3) {
      ii++;
    }
  }

  emxInit_int32_T(&r0, 1);
  i1 = b_ii->size[0];
  if (1 > idx) {
    b_ii->size[0] = 0;
  } else {
    b_ii->size[0] = idx;
  }

  emxEnsureCapacity((emxArray__common *)b_ii, i1, (int)sizeof(int));
  i1 = r0->size[0];
  r0->size[0] = 1 + b_ii->size[0];
  emxEnsureCapacity((emxArray__common *)r0, i1, (int)sizeof(int));
  r0->data[0] = 1;
  ii = b_ii->size[0];
  for (i1 = 0; i1 < ii; i1++) {
    r0->data[i1 + 1] = b_ii->data[i1] + 1;
  }

  emxInit_real_T1(&iTemp, 1);
  i1 = iTemp->size[0];
  iTemp->size[0] = r0->size[0];
  emxEnsureCapacity((emxArray__common *)iTemp, i1, (int)sizeof(double));
  ii = r0->size[0];
  for (i1 = 0; i1 < ii; i1++) {
    iTemp->data[i1] = 1.0 + (double)(r0->data[i1] - 1);
  }

  emxFree_int32_T(&r0);
  emxInit_real_T1(&c_yTemp, 1);
  i1 = c_yTemp->size[0];
  c_yTemp->size[0] = iTemp->size[0];
  emxEnsureCapacity((emxArray__common *)c_yTemp, i1, (int)sizeof(double));
  ii = iTemp->size[0];
  for (i1 = 0; i1 < ii; i1++) {
    c_yTemp->data[i1] = b_yTemp[(int)iTemp->data[i1] - 1];
  }

  emxInit_real_T1(&s, 1);
  emxInit_boolean_T(&b_x, 1);
  emxInit_real_T1(&r1, 1);
  b_diff(c_yTemp, s);
  b_sign(s);
  b_diff(s, r1);
  i1 = b_x->size[0];
  b_x->size[0] = r1->size[0];
  emxEnsureCapacity((emxArray__common *)b_x, i1, (int)sizeof(boolean_T));
  ii = r1->size[0];
  emxFree_real_T(&c_yTemp);
  for (i1 = 0; i1 < ii; i1++) {
    b_x->data[i1] = (r1->data[i1] < 0.0);
  }

  emxFree_real_T(&r1);
  nx = b_x->size[0];
  idx = 0;
  i1 = b_ii->size[0];
  b_ii->size[0] = b_x->size[0];
  emxEnsureCapacity((emxArray__common *)b_ii, i1, (int)sizeof(int));
  ii = 1;
  exitg2 = false;
  while ((!exitg2) && (ii <= nx)) {
    guard2 = false;
    if (b_x->data[ii - 1]) {
      idx++;
      b_ii->data[idx - 1] = ii;
      if (idx >= nx) {
        exitg2 = true;
      } else {
        guard2 = true;
      }
    } else {
      guard2 = true;
    }

    if (guard2) {
      ii++;
    }
  }

  if (b_x->size[0] == 1) {
    if (idx == 0) {
      i1 = b_ii->size[0];
      b_ii->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)b_ii, i1, (int)sizeof(int));
    }
  } else {
    i1 = b_ii->size[0];
    if (1 > idx) {
      b_ii->size[0] = 0;
    } else {
      b_ii->size[0] = idx;
    }

    emxEnsureCapacity((emxArray__common *)b_ii, i1, (int)sizeof(int));
  }

  if (1.0 > (double)s->size[0] - 1.0) {
    ii = 0;
  } else {
    ii = (int)((double)s->size[0] - 1.0);
  }

  if (2 > s->size[0]) {
    i1 = 0;
  } else {
    i1 = 1;
  }

  idx = b_x->size[0];
  b_x->size[0] = ii;
  emxEnsureCapacity((emxArray__common *)b_x, idx, (int)sizeof(boolean_T));
  for (idx = 0; idx < ii; idx++) {
    b_x->data[idx] = (s->data[idx] != s->data[i1 + idx]);
  }

  emxFree_real_T(&s);
  emxInit_int32_T(&c_ii, 1);
  nx = b_x->size[0];
  idx = 0;
  i1 = c_ii->size[0];
  c_ii->size[0] = b_x->size[0];
  emxEnsureCapacity((emxArray__common *)c_ii, i1, (int)sizeof(int));
  ii = 1;
  exitg1 = false;
  while ((!exitg1) && (ii <= nx)) {
    guard1 = false;
    if (b_x->data[ii - 1]) {
      idx++;
      c_ii->data[idx - 1] = ii;
      if (idx >= nx) {
        exitg1 = true;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1) {
      ii++;
    }
  }

  if (b_x->size[0] == 1) {
    if (idx == 0) {
      i1 = c_ii->size[0];
      c_ii->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)c_ii, i1, (int)sizeof(int));
    }
  } else {
    i1 = c_ii->size[0];
    if (1 > idx) {
      c_ii->size[0] = 0;
    } else {
      c_ii->size[0] = idx;
    }

    emxEnsureCapacity((emxArray__common *)c_ii, i1, (int)sizeof(int));
  }

  emxFree_boolean_T(&b_x);
  i1 = iInflect->size[0];
  iInflect->size[0] = c_ii->size[0];
  emxEnsureCapacity((emxArray__common *)iInflect, i1, (int)sizeof(double));
  ii = c_ii->size[0];
  for (i1 = 0; i1 < ii; i1++) {
    iInflect->data[i1] = iTemp->data[c_ii->data[i1]] - 1.0;
  }

  emxFree_int32_T(&c_ii);
  i1 = iPk->size[0];
  iPk->size[0] = b_ii->size[0];
  emxEnsureCapacity((emxArray__common *)iPk, i1, (int)sizeof(double));
  ii = b_ii->size[0];
  for (i1 = 0; i1 < ii; i1++) {
    iPk->data[i1] = iTemp->data[b_ii->data[i1]] - 1.0;
  }

  emxFree_int32_T(&b_ii);
  emxFree_real_T(&iTemp);
}

//
// Arguments    : const double Yin[49]
//                emxArray_real_T *Ypk
//                emxArray_real_T *Xpk
// Return Type  : void
//
static void findpeaks(const double Yin[49], emxArray_real_T *Ypk,
                      emxArray_real_T *Xpk)
{
  boolean_T x[49];
  int cdiff;
  emxArray_int32_T *ii;
  int idx;
  int k;
  boolean_T exitg1;
  emxArray_real_T *iInfite;
  boolean_T guard1 = false;
  double yTemp[49];
  emxArray_real_T *iPk;
  emxArray_real_T *iInflect;
  int ndbl;
  emxArray_real_T *base;
  double b_idx;
  int apnd;
  emxArray_real_T *y;
  emxArray_real_T *c_idx;
  emxArray_real_T *d_idx;
  for (cdiff = 0; cdiff < 49; cdiff++) {
    x[cdiff] = (rtIsInf(Yin[cdiff]) && (Yin[cdiff] > 0.0));
  }

  emxInit_int32_T(&ii, 1);
  idx = 0;
  k = ii->size[0];
  ii->size[0] = 49;
  emxEnsureCapacity((emxArray__common *)ii, k, (int)sizeof(int));
  cdiff = 1;
  exitg1 = false;
  while ((!exitg1) && (cdiff < 50)) {
    guard1 = false;
    if (x[cdiff - 1]) {
      idx++;
      ii->data[idx - 1] = cdiff;
      if (idx >= 49) {
        exitg1 = true;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1) {
      cdiff++;
    }
  }

  emxInit_real_T1(&iInfite, 1);
  k = ii->size[0];
  if (1 > idx) {
    ii->size[0] = 0;
  } else {
    ii->size[0] = idx;
  }

  emxEnsureCapacity((emxArray__common *)ii, k, (int)sizeof(int));
  k = iInfite->size[0];
  iInfite->size[0] = ii->size[0];
  emxEnsureCapacity((emxArray__common *)iInfite, k, (int)sizeof(double));
  cdiff = ii->size[0];
  for (k = 0; k < cdiff; k++) {
    iInfite->data[k] = ii->data[k];
  }

  memcpy(&yTemp[0], &Yin[0], 49U * sizeof(double));
  k = ii->size[0];
  ii->size[0] = iInfite->size[0];
  emxEnsureCapacity((emxArray__common *)ii, k, (int)sizeof(int));
  cdiff = iInfite->size[0];
  for (k = 0; k < cdiff; k++) {
    ii->data[k] = (int)iInfite->data[k];
  }

  cdiff = ii->size[0];
  for (k = 0; k < cdiff; k++) {
    yTemp[ii->data[k] - 1] = rtNaN;
  }

  emxInit_real_T1(&iPk, 1);
  emxInit_real_T1(&iInflect, 1);
  findLocalMaxima(yTemp, iPk, iInflect);
  if (!(iPk->size[0] == 0)) {
    idx = iPk->size[0] - 1;
    ndbl = 0;
    for (cdiff = 0; cdiff <= idx; cdiff++) {
      if (Yin[(int)iPk->data[cdiff] - 1] > rtMinusInf) {
        ndbl++;
      }
    }

    k = 0;
    for (cdiff = 0; cdiff <= idx; cdiff++) {
      if (Yin[(int)iPk->data[cdiff] - 1] > rtMinusInf) {
        iPk->data[k] = iPk->data[cdiff];
        k++;
      }
    }

    k = iPk->size[0];
    iPk->size[0] = ndbl;
    emxEnsureCapacity((emxArray__common *)iPk, k, (int)sizeof(double));
  }

  cdiff = iPk->size[0];
  emxInit_real_T1(&base, 1);
  k = base->size[0];
  base->size[0] = cdiff;
  emxEnsureCapacity((emxArray__common *)base, k, (int)sizeof(double));
  for (k = 0; k + 1 <= cdiff; k++) {
    if ((Yin[(int)(iPk->data[k] - 1.0) - 1] >= Yin[(int)(iPk->data[k] + 1.0) - 1])
        || rtIsNaN(Yin[(int)(iPk->data[k] + 1.0) - 1])) {
      b_idx = Yin[(int)(iPk->data[k] - 1.0) - 1];
    } else {
      b_idx = Yin[(int)(iPk->data[k] + 1.0) - 1];
    }

    base->data[k] = b_idx;
  }

  idx = iPk->size[0] - 1;
  ndbl = 0;
  for (cdiff = 0; cdiff <= idx; cdiff++) {
    if (Yin[(int)iPk->data[cdiff] - 1] - base->data[cdiff] >= 0.0) {
      ndbl++;
    }
  }

  k = 0;
  for (cdiff = 0; cdiff <= idx; cdiff++) {
    if (Yin[(int)iPk->data[cdiff] - 1] - base->data[cdiff] >= 0.0) {
      iPk->data[k] = iPk->data[cdiff];
      k++;
    }
  }

  k = iPk->size[0];
  iPk->size[0] = ndbl;
  emxEnsureCapacity((emxArray__common *)iPk, k, (int)sizeof(double));
  combinePeaks(iPk, iInfite, base);
  emxFree_real_T(&iInfite);
  if (base->size[0] < 1) {
    ndbl = 0;
    apnd = 0;
  } else {
    ndbl = (int)std::floor(((double)base->size[0] - 1.0) + 0.5);
    apnd = ndbl + 1;
    cdiff = (ndbl - base->size[0]) + 1;
    idx = base->size[0];
    if (std::abs((double)cdiff) < 4.4408920985006262E-16 * (double)idx) {
      ndbl++;
      apnd = base->size[0];
    } else if (cdiff > 0) {
      apnd = ndbl;
    } else {
      ndbl++;
    }
  }

  emxInit_real_T(&y, 2);
  k = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = ndbl;
  emxEnsureCapacity((emxArray__common *)y, k, (int)sizeof(double));
  if (ndbl > 0) {
    y->data[0] = 1.0;
    if (ndbl > 1) {
      y->data[ndbl - 1] = apnd;
      idx = (ndbl - 1) / 2;
      for (k = 1; k < idx; k++) {
        y->data[k] = 1.0 + (double)k;
        y->data[(ndbl - k) - 1] = apnd - k;
      }

      if (idx << 1 == ndbl - 1) {
        y->data[idx] = (1.0 + (double)apnd) / 2.0;
      } else {
        y->data[idx] = 1.0 + (double)idx;
        y->data[idx + 1] = apnd - idx;
      }
    }
  }

  emxInit_real_T1(&c_idx, 1);
  k = c_idx->size[0];
  c_idx->size[0] = y->size[1];
  emxEnsureCapacity((emxArray__common *)c_idx, k, (int)sizeof(double));
  cdiff = y->size[1];
  for (k = 0; k < cdiff; k++) {
    c_idx->data[k] = y->data[y->size[0] * k];
  }

  emxFree_real_T(&y);
  if (c_idx->size[0] == 0) {
  } else {
    k = iInflect->size[0];
    iInflect->size[0] = c_idx->size[0];
    emxEnsureCapacity((emxArray__common *)iInflect, k, (int)sizeof(double));
    cdiff = c_idx->size[0];
    for (k = 0; k < cdiff; k++) {
      iInflect->data[k] = Yin[(int)base->data[(int)c_idx->data[k] - 1] - 1];
    }

    emxInit_real_T1(&d_idx, 1);
    sort(iInflect, ii);
    k = d_idx->size[0];
    d_idx->size[0] = ii->size[0];
    emxEnsureCapacity((emxArray__common *)d_idx, k, (int)sizeof(double));
    cdiff = ii->size[0];
    for (k = 0; k < cdiff; k++) {
      d_idx->data[k] = c_idx->data[ii->data[k] - 1];
    }

    k = c_idx->size[0];
    c_idx->size[0] = d_idx->size[0];
    emxEnsureCapacity((emxArray__common *)c_idx, k, (int)sizeof(double));
    cdiff = d_idx->size[0];
    for (k = 0; k < cdiff; k++) {
      c_idx->data[k] = d_idx->data[k];
    }

    emxFree_real_T(&d_idx);
  }

  emxFree_int32_T(&ii);
  emxFree_real_T(&iInflect);
  if (c_idx->size[0] > 1) {
    b_idx = c_idx->data[0];
    k = c_idx->size[0];
    c_idx->size[0] = 1;
    emxEnsureCapacity((emxArray__common *)c_idx, k, (int)sizeof(double));
    c_idx->data[0] = b_idx;
  }

  k = iPk->size[0];
  iPk->size[0] = c_idx->size[0];
  emxEnsureCapacity((emxArray__common *)iPk, k, (int)sizeof(double));
  cdiff = c_idx->size[0];
  for (k = 0; k < cdiff; k++) {
    iPk->data[k] = base->data[(int)c_idx->data[k] - 1];
  }

  emxFree_real_T(&base);
  emxFree_real_T(&c_idx);
  k = Ypk->size[0];
  Ypk->size[0] = iPk->size[0];
  emxEnsureCapacity((emxArray__common *)Ypk, k, (int)sizeof(double));
  cdiff = iPk->size[0];
  for (k = 0; k < cdiff; k++) {
    Ypk->data[k] = Yin[(int)iPk->data[k] - 1];
  }

  k = Xpk->size[0];
  Xpk->size[0] = iPk->size[0];
  emxEnsureCapacity((emxArray__common *)Xpk, k, (int)sizeof(double));
  cdiff = iPk->size[0];
  for (k = 0; k < cdiff; k++) {
    Xpk->data[k] = 1.0 + (double)((int)iPk->data[k] - 1);
  }

  emxFree_real_T(&iPk);
}

//
// Arguments    : const double y[50]
//                emxArray_real_T *iPk
//                emxArray_real_T *iInf
//                emxArray_real_T *iInflect
// Return Type  : void
//
static void getAllPeaks(const double y[50], emxArray_real_T *iPk,
  emxArray_real_T *iInf, emxArray_real_T *iInflect)
{
  boolean_T x[50];
  int ii;
  emxArray_int32_T *b_ii;
  int idx;
  int i3;
  boolean_T exitg1;
  boolean_T guard1 = false;
  double yTemp[50];
  for (ii = 0; ii < 50; ii++) {
    x[ii] = (rtIsInf(y[ii]) && (y[ii] > 0.0));
  }

  emxInit_int32_T(&b_ii, 1);
  idx = 0;
  i3 = b_ii->size[0];
  b_ii->size[0] = 50;
  emxEnsureCapacity((emxArray__common *)b_ii, i3, (int)sizeof(int));
  ii = 1;
  exitg1 = false;
  while ((!exitg1) && (ii < 51)) {
    guard1 = false;
    if (x[ii - 1]) {
      idx++;
      b_ii->data[idx - 1] = ii;
      if (idx >= 50) {
        exitg1 = true;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1) {
      ii++;
    }
  }

  i3 = b_ii->size[0];
  if (1 > idx) {
    b_ii->size[0] = 0;
  } else {
    b_ii->size[0] = idx;
  }

  emxEnsureCapacity((emxArray__common *)b_ii, i3, (int)sizeof(int));
  i3 = iInf->size[0];
  iInf->size[0] = b_ii->size[0];
  emxEnsureCapacity((emxArray__common *)iInf, i3, (int)sizeof(double));
  ii = b_ii->size[0];
  for (i3 = 0; i3 < ii; i3++) {
    iInf->data[i3] = b_ii->data[i3];
  }

  memcpy(&yTemp[0], &y[0], 50U * sizeof(double));
  i3 = b_ii->size[0];
  b_ii->size[0] = iInf->size[0];
  emxEnsureCapacity((emxArray__common *)b_ii, i3, (int)sizeof(int));
  ii = iInf->size[0];
  for (i3 = 0; i3 < ii; i3++) {
    b_ii->data[i3] = (int)iInf->data[i3];
  }

  ii = b_ii->size[0];
  for (i3 = 0; i3 < ii; i3++) {
    yTemp[b_ii->data[i3] - 1] = rtNaN;
  }

  emxFree_int32_T(&b_ii);
  b_findLocalMaxima(yTemp, iPk, iInflect);
}

//
// Arguments    : emxArray_real_T *idx
// Return Type  : void
//
static void keepAtMostNpPeaks(emxArray_real_T *idx)
{
  emxArray_real_T *b_idx;
  int i10;
  int loop_ub;
  if (idx->size[0] > 50) {
    emxInit_real_T1(&b_idx, 1);
    i10 = b_idx->size[0];
    b_idx->size[0] = 50;
    emxEnsureCapacity((emxArray__common *)b_idx, i10, (int)sizeof(double));
    for (i10 = 0; i10 < 50; i10++) {
      b_idx->data[i10] = idx->data[i10];
    }

    i10 = idx->size[0];
    idx->size[0] = b_idx->size[0];
    emxEnsureCapacity((emxArray__common *)idx, i10, (int)sizeof(double));
    loop_ub = b_idx->size[0];
    for (i10 = 0; i10 < loop_ub; i10++) {
      idx->data[i10] = b_idx->data[i10];
    }

    emxFree_real_T(&b_idx);
  }
}

//
// KNNCLASSIFY - Simple K Nearest Neighbor Algorithm using 2norm method:
//  Classify using the Nearest neighbor algorithm
//  Inputs:
//   tX          - Train samples
//  tY          - Train labels
//    testSamples - Test  samples
//  KNN      - Number of nearest neighbors
//
//  Outputs
//  Y           - Predicted targets
// Arguments    : const emxArray_real_T *testSamples
// Return Type  : double
//
static double knnclassify(const emxArray_real_T *testSamples)
{
  double Y;
  emxArray_real_T *Uc;
  int low_ip1;
  int k;
  int nb;
  static const signed char a[245] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 };

  emxArray_real_T *b_testSamples;
  int nbins;
  int exitg5;
  int i7;
  boolean_T eok;
  emxArray_real_T *r4;
  static double x[25480];
  static double y[25480];
  int b_k;
  static const double tX[25480] = { 8.09154631093208E-5, 0.000244302455926219,
    0.104682824331086, 0.324052425153446, 0.0160228177622438, 0.106634131845299,
    0.055798992146868, 0.0330228453481924, 0.0233441111070391,
    0.00504321184494825, 0.01065438799826, 0.00585703467583621,
    0.0344886600837497, 0.00772742672694014, 0.0567746459039746,
    0.0315601427458316, 0.0440040073320862, 0.026028325989031,
    0.000163386992816898, 8.09154631093208E-5, 8.09154631093208E-5,
    0.000163386992816898, 0.000161830926218642, 0.000242746389327963,
    0.12038353630749, 0.264674479830588, 0.0300958840768726, 0.0866246714583227,
    0.0160243738288421, 0.0273291976651731, 0.0570189483599008,
    0.0601886560205487, 0.00154673019866663, 0.00170856112488527,
    0.000650435838071076, 0.00113748468332527, 0.134942095400775,
    0.134371018959215, 0.204241521354114, 0.00984212123397028,
    0.000812266764289721, 0.00121995621303284, 0.00187039205110392,
    0.0017894765879946, 0.000487048845254183, 0.00097565375710662, 0.0,
    0.000325217919035538, 0.00121995621303284, 0.0130958564909239,
    0.0194399400120143, 0.000407689448743111, 0.00309034826413675,
    0.0181375122692739, 0.00577611921272691, 0.000244302455926215,
    0.0047179939259127, 0.0116315978219649, 0.0305020174590175,
    0.00764651126383081, 0.051000082757846, 0.00252082788917501,
    0.00740220880790461, 0.00243835635946742, 0.0258649389962141,
    0.00154673019866662, 0.00341556618317229, 0.000490160978450699,
    0.00048704884525419, 0.0581579891098244, 0.00113748468332525,
    0.0388814360906269, 0.000813822830887989, 0.00138334320584974,
    0.00406755808784162, 0.0204171498357192, 0.0627935115060295,
    0.0153723819241727, 0.00292851733791811, 0.00317126372724607,
    0.0074022088079046, 0.253126909604929, 0.0701957203139341, 0.160480260411355,
    0.0155373249835877, 0.000813822830888045, 0.00634563958768863,
    0.00910921386619168, 0.0184642862549076, 0.0627141521095185,
    0.16153838569817, 0.369439775691381, 0.0624682935869938, 0.0470959116628211,
    0.164222600580161, 0.115419683859053, 0.0109811619838938, 0.0319662761279764,
    0.196679037686588, 0.148605916200063, 0.192284705613113, 0.0117125132850742,
    0.21082990733113, 0.228888060203893, 0.640785113028694, 0.00634408352109039,
    0.181547846085146, 0.0639325522559531, 0.120297952644586, 0.311852863023117,
    0.0990716481777738, 0.0827220564298963, 0.108098390514258,
    0.00341712224977053, 0.27354250337405, 0.0165923381372055,
    0.0732035970483632, 0.328446757226921, 0.211074209787057, 0.195457525406957,
    0.555869002695257, 0.0707652406888959, 0.0497801265448132, 0.281839450475952,
    0.317221292787101, 0.148604360133464, 0.353336042466028, 0.103951473029905,
    0.188624836974015, 0.229619411505073, 0.0634439473441006, 0.0353833983777471,
    0.275008318109608, 0.303554359854617, 0.0641768547118793, 0.012932469498107,
    0.177642118923523, 0.290136397577854, 0.0197667139976481, 0.0141539817777381,
    0.150314477324948, 0.334302235836159, 0.294526061451535, 0.273786805829977,
    0.272320991094419, 0.0966301796851099, 0.0273291976651731,
    0.0954102234720771, 0.229375109049147, 0.207170038692032, 0.0217164654452632,
    0.0658854158367643, 0.119079552498152, 0.133966441643668, 0.0551470002421987,
    0.0287950124007305, 0.535859542308281, 0.0566112589111577, 0.196431623097466,
    0.0546599513969445, 0.165199810403866, 0.0278178025770257,
    0.0712522895341501, 0.108098390514258, 0.278909377071436, 0.0307463199149437,
    0.0551501123753951, 0.280862240652247, 0.0893119984735112, 0.204241521354114,
    0.346501797966487, 0.167151117918079, 0.0812577977609372, 0.10858699542611,
    0.219859761800811, 0.102733072883471, 0.00780834219004956, 0.169102425432293,
    0.254021647898927, 0.0541713464850921, 0.385302318594005, 0.0646639035571334,
    0.40628587667149, 0.288916441364821, 0.115663986314979, 0.00902829840308228,
    0.144946047560964, 0.100780209302659, 1.37381229793542, 0.250604525649156,
    1.54755024576391, 1.75691745049268, 0.81379326562261, 0.722287213245558,
    1.11564217648554, 0.610771700481529, 1.27644921088253, 1.35282562772474,
    0.373589805308931, 0.985338715614165, 0.447282007269147, 1.23252611901355,
    0.056611258911158, 0.450454827062991, 0.696422274249343, 0.977528817357517,
    0.826238686275463, 0.793295200323782, 0.794272410147486, 1.4260307808397,
    0.129816412026119, 0.664943046966621, 1.22715768924957, 1.4499444123217,
    1.92675278329275, 0.848443756632578, 0.806229225888486, 0.675435604038662,
    0.641763878918997, 1.39260024604276, 1.2127594050159, 1.41016979400367,
    1.51143705215158, 0.463874345406352, 1.38869607494774, 0.628585550898366,
    2.12489296751511, 0.427516849338096, 0.696422274249344, 1.36429539462048,
    0.621754418532022, 0.0797935390919786, 0.40409026670135, 0.491936450439304,
    1.21715218102278, 0.179840841026858, 1.09660836985567, 0.700082142888442,
    0.567092911068479, 1.40357985196006, 0.924821729541384, 0.793296756390379,
    0.552942041423937, 0.668360169216391, 0.363584297082144,
    -1.01144328886651E-5, -0.000244302455926219, -0.0523414121655428,
    -0.0810131062883614, -0.0160228177622438, -0.0266585329613247,
    0.055798992146868, 0.0330228453481924, -0.0233441111070391,
    0.00252160592247412, -0.000591910444347777, 0.000167343847881035,
    -0.0114962200279166, -0.00045455451334942, -0.00236561024599894,
    0.0315601427458316, -0.0440040073320862, -0.00520566519780619,
    0.000163386992816898, 8.09154631093208E-5, 2.31187037455202E-6,
    -0.000163386992816898, 0.000161830926218642, -3.46780556182804E-5,
    -0.0601917681537452, -0.264674479830588, -0.0300958840768726,
    -0.00666343626602483, 0.000572299065315788, -0.0136645988325866,
    0.00633543870665565, 0.0150471640051372, -7.73365099333318E-5,
    0.00170856112488527, -0.000162608959517769, -0.000103407698484115,
    0.0674710477003875, -0.134371018959215, 0.204241521354114,
    0.00164035353899505, 0.000812266764289721, -0.00121995621303284,
    -8.90662881478056E-5, 6.62769106664665E-5, -2.21385838751901E-5,
    -6.96895540790443E-5, 0.0, 0.000325217919035538, -2.9046576500782E-5,
    -0.0130958564909239, -0.00971997000600717, -0.000101922362185778,
    0.000772587066034189, -0.0181375122692739, 0.00288805960636345,
    -0.000244302455926215, -0.0047179939259127, -0.0116315978219649,
    -0.00762550436475436, -0.00764651126383081, 0.025500041378923,
    0.000148283993480883, 0.000616850733992051, -0.000270928484385269,
    0.012932469498107, 0.000103115346577775, 0.000170778309158615,
    1.48533629833545E-5, -0.00048704884525419, -0.0290789945549122,
    0.000568742341662626, -0.00388814360906269, 0.000406911415443995,
    -0.00138334320584974, 0.000813511617568324, -0.0204171498357192,
    0.0627935115060295, 0.000614895276966908, -0.000585703467583623,
    -0.00158563186362303, -0.0012337014679841, -0.0843756365349765,
    -0.0701957203139341, 0.160480260411355, -0.00776866249179387,
    -0.000813822830888045, 0.000192292108717837, -0.00910921386619168,
    -0.0184642862549076, 0.00313570760547592, 0.0807691928490848,
    -0.0335854341537619, 0.00416455290579959, -0.00428144651480192,
    -0.0547408668600538, 0.0115419683859053, -0.000915096831991147,
    0.00399578451599705, 0.0245848797108236, -0.148605916200063,
    -0.00468987086861252, 0.00146406416063427, -0.21082990733113,
    0.00635800167233036, -0.0194177306978392, -0.00634408352109039,
    0.181547846085146, 0.00376073836799724, 0.0240595905289173,
    -0.00842845575738155, 0.0123839560222217, -0.0275740188099654,
    0.00251391605847111, -0.000113904074992351, -0.0248675003067318,
    0.000638166851430981, -0.00665487245894211, 0.0117302413295329,
    -0.0175895174822547, 0.00698062590739134, 0.0794098575278938,
    0.0101093200984137, 0.00121414942792227, -0.00909159517664362,
    -0.317221292787101, -0.00675474364243019, -0.353336042466028,
    0.00358453355275535, 0.0188624836974015, -0.229619411505073,
    -0.00288381578836821, -0.00353833983777471, 0.137504159054804,
    0.0433649085506596, -0.00534807122598994, -0.000417176435422807,
    0.0253774455605032, -0.290136397577854, 0.00073210051843141,
    -0.00117949848147817, -0.00715783225356893, 0.0167151117918079,
    -0.0267750964955941, -0.00912622686099922, -0.0247564537358563,
    -0.00603938623031937, -0.00118822598544231, 0.00238525558680193,
    -0.0191145924207623, 0.0295957198131474, -0.000723882181508773,
    -0.00549045131973036, 0.00661553069434176, -0.0446554805478894,
    0.0183823334140662, -0.00179968827504565, -0.535859542308281,
    -0.0188704196370526, 0.0122769764435916, 0.00780856448527778,
    -0.165199810403866, -0.0278178025770257, -0.00203577970097572,
    -0.108098390514258, -0.0232424480892863, 0.0102487733049812,
    0.0183833707917984, -0.0351077800815309, -0.00307972408529349,
    0.0136161014236076, 0.0138600719186595, -0.167151117918079,
    0.011608256822991, 0.0361956651420368, 0.0104695124667053,
    0.00428054470347794, 0.000251882006130631, -0.0153729477665721,
    0.0362888068427038, 0.00773876378358459, -0.0132862868480691,
    0.00923770050816191, -0.40628587667149, -0.0262651310331656,
    -0.115663986314979, -0.000531076376651899, -0.144946047560964,
    0.0503901046513295, 1.37381229793542, 0.250604525649156, -0.140686385978538,
    -0.0488032625136856, 0.058128090401615, 0.722287213245558, 1.11564217648554,
    0.610771700481529, 1.27644921088253, 1.35282562772474, 0.026684986093495,
    0.0447881234370075, -0.0372735006057622, 1.23252611901355,
    -0.00314506993950878, 0.450454827062991, 0.0497444481606674,
    0.977528817357517, -0.0217431233230385, 0.793295200323782, 0.794272410147486,
    0.0365648918164025, 0.00927260085900853, -0.0604493679060565,
    -0.111559789931779, 1.4499444123217, 1.92675278329275, 0.848443756632578,
    0.0310088163803264, -0.0182550163253693, -0.0267401616216249,
    1.39260024604276, 0.044917015000589, 0.054237299769372, 1.51143705215158,
    0.463874345406352, -0.126245097722522, 0.628585550898366, -0.193172087955919,
    0.427516849338096, 0.696422274249344, 1.36429539462048, -0.036573789325413,
    0.0797935390919786, 0.40409026670135, 0.061492056304913, -0.0553250991373991,
    0.179840841026858, -0.0645063746973923, 0.700082142888442,
    -0.0257769505031127, 1.40357985196006, 0.071140133041645, -0.07211788694458,
    0.552942041423937, 0.668360169216391, 0.030298691423512, 0.711147332468642,
    0.703012216292959, 0.275663422147473, 3.86610814996478, 2.17743822440502,
    -2.83620524304316, -3.6731745645802, -1.54608443102836, -1.89160590703794,
    -3.48527952284077, -5.1966477918693, -4.77327941608215, -3.55124741020724,
    -2.54395881915807, -2.59632668445579, -4.34641922684852, -2.99684733438074,
    -1.49321395621941, 0.32689535882846, 0.325997508401266, 0.327710737725946,
    0.325760986278331, 0.325599155352112, 0.332673034107785, 0.348124775428469,
    -1.8785224990798, -1.70933604605121, -7.16171871514411, -9.37298738346182,
    -3.53058907004879, -5.27814989208616, -4.92798511334738, -4.95686882154421,
    -5.49280149898261, -5.78863154424348, -5.95651556952934, -6.38866482913044,
    -7.29567426432159, -6.53784337783867, -1.95213067344371, -2.40917388648372,
    -2.5174356639908, -2.48660220434635, -2.58811998921658, -2.62634943340254,
    -2.61944049770628, -2.63252701779762, -2.70426324404383, -2.7336402253523,
    -2.81611331112648, -1.59034052114936, -1.57154790484222, -1.75505017057478,
    -2.42235221450435, -4.156812511851, -4.34528796643159, -4.32225973684399,
    -4.37513176771954, -5.86045957841899, -6.64618762687545, -6.48813949855717,
    -5.86346745515341, -5.74251595453756, -3.2925949081119, -3.40401861294663,
    -3.43680804830509, -3.65658378250959, -3.82886058380233, -3.92369506263305,
    -3.72808348663287, -3.86748838103744, -3.91337678502001, -4.1117566034985,
    -4.30224872039042, -5.05129561092619, -1.53242061022907, -0.450856292245341,
    -0.829335258805992, -0.525611287692165, -0.35813495185505,
    -0.393186908047369, 2.93795643790313, 41.8963555824798, 37.475926716085,
    37.5004845591387, 34.9770378098788, 30.8374649598699, 30.0549003948077,
    31.4333497712059, -10.7879529747213, -9.77456149047379, -10.2833299133064,
    -11.6917942583184, -10.7789262323848, -12.8269588463457, -10.4658518570821,
    -11.3735926435744, -11.1691068197644, -12.7015258739268, -10.0132045321821,
    -11.5134207880937, -10.6539880891442, -11.4655857447967, -10.8113826694913,
    -11.4582691196517, -11.6905774142386, -12.4514177335225, -10.4087551053923,
    -13.0009348723637, -13.3279220391215, -12.0885694559414, -11.6964251125148,
    -9.915354396284, -11.1866732555921, -9.63399577038691, -10.5715150033701,
    -12.0343903291233, -10.7135340896597, -11.9365401932251, -10.4980655477992,
    -11.9667963521616, -9.79285927760268, -10.5222234817371, -11.4834073755465,
    -10.8133339770055, -11.7344942818412, -10.3941109626361, -11.2418280361673,
    -11.9943682962161, -11.7115531919831, -11.4887711371107, -10.6630257239469,
    -10.8291965199081, -9.75674452792375, -11.4150789351505, -10.7379316578538,
    -11.9111638591408, -10.4397488398963, -10.851398478132, -11.4358197468387,
    -11.6876457847675, -10.9975738183056, -11.0095275219134, -11.1974163393864,
    -10.1437554076426, -11.8882274374825, -11.2096127893836, -11.075901542662,
    -12.1747024103546, -9.84776197538895, -11.1825278941743, -12.0822160360207,
    -11.1608145408623, -10.0478550231921, -12.1995870273939, -11.8325840519954,
    -10.3104101403159, -12.1947056464752, -10.8045499810583, -10.3811753810048,
    -11.1818027671395, -11.6080934359982, -10.9302179195335, -12.30696184694,
    -12.3562471443066, -11.2162042874938, -12.5497549182659, -11.2720748587042,
    -11.6869097652665, -12.1742122493762, -11.2547589495988, -11.3518761781291,
    -11.1156652685139, -12.0592796143624, -11.8240474706374, -11.5221923355081,
    -12.0299975531164, -10.9507159848323, -11.5680729591577, -11.0236815036911,
    -12.1390731534544, -11.6400597121262, -11.5412354785375, -11.2345051867559,
    -11.1224995130134, -11.7925651312215, -13.1434410074986, -11.1917973829001,
    -4.44353489931228, -9.23186147939957, -10.1491207252734, -14.3008091098164,
    -11.8789501684237, -13.375259141107, -10.4529287239836, -14.8759624458638,
    -17.8856582221464, -13.4570055437798, -16.8764058750502, -14.6243775982578,
    -10.8086984546093, -12.5385294538261, -14.8876671788159, -12.7020238152383,
    -16.3393170401293, -15.164877331162, -18.3875908444129, -15.6897525993532,
    -14.7327171790948, -18.8670896545323, -13.0834157384709, -13.1717520831873,
    -17.318313228289, -17.525483266981, -17.5137754218957, -15.1558474766924,
    -14.4181816572562, -15.8017598291623, -15.73099614454, -15.0955747930755,
    -13.7942347408873, -15.9845229632607, -18.9283395479729, -15.0257836500771,
    -17.6655494896898, -18.4183371643278, -17.6323710376818, -20.0308033964378,
    -23.3052561240806, -22.3787242773477, -21.129605820197, -19.5991302932157,
    -20.5493345811079, -18.0769641618692, -15.0948449978409, -26.4606292992253,
    -18.1562659839161, -20.5949662341017, -24.0153344387963, -26.4533080058805,
    -20.3892604540787, -18.3968696695383, -24.7386144225315, -17.6067441768751,
    -8.55886734385357, 0.000457119255178844, 0.000455345452249949,
    0.0661455085851307, 0.176083329377243, 0.07609179831709, 0.174000277129822,
    0.0425038680588661, 0.0552608186077743, 0.0136653192637311,
    0.0151812903535802, 0.019951886183515, 0.0171155660254629, 0.041049116066785,
    0.0134495587228131, 0.0295507403963712, 0.0241299651789822,
    0.0452643510039559, 0.064215973161458, 0.000346690224156388,
    0.000325480372213072, 0.000298602149516788, 0.000257992836793177,
    0.000299029678758094, 0.000347696697164504, 0.0278522024416967,
    0.0713602088312158, 0.0454824206336368, 0.0748257791116681,
    0.0754969731074406, 0.0229184436934307, 0.127643953878272,
    0.0555566225372506, 0.00910446542038037, 0.00642266450795653,
    0.00302152652512267, 0.00273080123594629, 0.0280752532376485,
    0.14500481488212, 0.188709579784811, 0.0129156435102175, 0.00459952256267937,
    0.00290672715246149, 0.00285572904603352, 0.00205510379866526,
    0.00208440150802401, 0.00243357382580788, 0.00289432190530794,
    0.00224756718497758, 0.00240995709890518, 0.00375339700978922,
    0.0102520636866222, 0.00320282519836017, 0.00292942713631234,
    0.0287867561001097, 0.0126585946416211, 0.00202972945920737,
    0.00390612795074279, 0.00690617551049398, 0.0282110139209453,
    0.00987196732197399, 0.0279694535575448, 0.0127832522563022,
    0.0174507775278575, 0.0154480424373721, 0.0176170719312207,
    0.00539939391739025, 0.00511515081321851, 0.00362476929073195,
    0.00386145343956773, 0.0168416506872694, 0.00469930968850507,
    0.00996946231513324, 0.00333397464256584, 0.00929167294660743,
    0.0254731501086212, 0.0309458814354578, 0.0381940968280991,
    0.00962526614873011, 0.00629771172959604, 0.0115190835776697,
    0.0103287328852102, 0.249785040331958, 0.109005378304134, 0.0711627149423633,
    0.0311029666479487, 0.0153633723535361, 0.0249387592430557,
    0.00851104266909615, 0.0829589843702074, 0.304303436187998, 0.27086217654548,
    0.30463714713226, 0.271414927980596, 0.285640367314156, 0.255853804276653,
    0.265460076059956, 0.24447622772188, 0.232852127696263, 0.242120085200691,
    0.225705508173295, 0.249933611337555, 0.249554500874019, 0.2732719458542,
    0.220413959786329, 0.254475412473186, 0.199294756698113, 0.207750222381285,
    0.20977704102324, 0.232117564794401, 0.230278106253545, 0.240971372839061,
    0.291541853033027, 0.247914883776663, 0.289755026592962, 0.299954101192666,
    0.281488847961773, 0.270411235951062, 0.274763627024325, 0.254739932176764,
    0.264621948517917, 0.279607437530276, 0.276103633718669, 0.274431854519229,
    0.271425133453621, 0.277114607767025, 0.289173938994287, 0.267261012878404,
    0.279881941809144, 0.249855836703524, 0.275546320398953, 0.237929600747118,
    0.290800293966661, 0.228462389724506, 0.269878301565214, 0.286459280488773,
    0.290412566604068, 0.329585465057627, 0.331570803421684, 0.343591105606845,
    0.303342021425385, 0.339064536421727, 0.344314763403221, 0.254576325018917,
    0.326330853661705, 0.295612180073202, 0.264690652593832, 0.293521417468903,
    0.339280929574816, 0.286893967616264, 0.311781183704759, 0.307526773251742,
    0.298607169558283, 0.249357735336094, 0.297635115129837, 0.264431064779331,
    0.316855923647177, 0.302195817895372, 0.289413982866116, 0.276494375358473,
    0.279630092655857, 0.281262730090356, 0.274212216780685, 0.28510451643012,
    0.27175836559797, 0.290847123573552, 0.290872496061649, 0.283798151698051,
    0.275565762928934, 0.280445854924815, 0.289221828501849, 0.331154125258437,
    0.290694623489805, 0.28285795747342, 0.291971237684621, 0.28555592889479,
    0.261009364820808, 0.267720102486467, 0.312520640785395, 0.292220461999606,
    0.293125124131938, 0.297900033574819, 0.299328963082431, 0.325943710989614,
    0.302041365491396, 0.283500201636736, 0.316475698436968, 0.271160954407959,
    0.306963228642845, 0.693044443069628, 0.712794896788199, 0.773239223580053,
    1.03341760033761, 0.796554053111978, 0.909603207691074, 0.649703013446223,
    0.904362335004985, 0.937298124129951, 0.970607653779548, 0.804598965089867,
    0.846346288451135, 0.757876035159117, 0.764595579776931, 0.896099765045261,
    0.870512139399764, 0.856361455954944, 0.835956322048036, 1.02188459547839,
    0.833352757582343, 0.70231613072067, 0.873223470429507, 0.918232235439818,
    0.775286590418254, 0.798278920744475, 0.876528122374015, 1.00274875178612,
    0.804765072742737, 0.816858937474218, 0.856764320206191, 0.844974288633314,
    0.972203714695043, 0.967725356384552, 0.765065719123715, 0.993128235704037,
    0.939532597645617, 1.12061735273473, 1.04896230391595, 0.927195796547378,
    1.03383268810218, 1.04436037168747, 0.982054652046149, 0.824853513477025,
    0.921506589815537, 1.0187078612906, 1.0616670771596, 1.03240447825867,
    1.01570079726729, 1.05769355327295, 0.888795723634297, 0.92625404046423,
    1.10668432655807, 1.02840961930002, 0.953635935992552, 1.3850438377011,
    1.31583071676714, 1.43239902958631, 0.0150471640051372, 0.0150471640051372,
    0.20440490834693, 0.465663821994347, 0.194480315583253, 0.502348092048236,
    -0.00666930144012594, 0.146816439612068, -0.0120377312041097,
    -0.0305829329221268, -0.0566112589111577, -0.0487220012579989,
    0.0113872953660387, -0.0116315978219649, 0.0133401589468501,
    0.0128515540349977, 0.116964857991122, 0.164060769653943,
    0.00732129334479528, 0.0072388218150877, 0.0072388218150877,
    0.00715790635197838, 0.0072388218150877, 0.0074022088079046,
    0.154218648419973, 0.257681516538025, 0.0994777815599187,
    0.00976120577086095, 0.0391241824799549, -0.0138272077921043,
    0.117534378366083, -0.0187895041739433, -0.0723104148209642,
    -0.0963049617660744, -0.110620774470031, -0.114688332557873,
    -0.022693675268968, 0.11208503313899, 0.250279307730121,
    -0.00422938901406027, -0.0375790083478865, -0.0447369146998649,
    -0.0427840511190534, -0.0463630042950426, -0.0480715654199279,
    -0.0456316529938622, -0.0475829605080754, -0.0498610420079225,
    -0.0492915216329607, -0.0454682660010453, 0.00032521791903554,
    -0.0237502444891839, -0.0293629767090939, 0.012932469498107,
    -0.0475829605080754, -0.0818273181358989, -0.0757259810041365,
    -0.0549042538528707, -0.0578327711907888, -0.0847542794072188,
    -0.060840647925218, -0.0934573598912656, -0.0617353862192153,
    -0.0181375122692739, -0.00845877802812054, -0.0562051255290129,
    -0.0593779453228572, -0.0659663312998738, -0.0676733363581608,
    -0.00911076993278987, -0.0667785980641635, -0.0425397486631272,
    -0.0714965919900762, -0.0734494555708877, -0.0563669564552315,
    0.083534323194186, 0.0376599238109958, 0.0090282984030823,
    0.00618225259487176, 0.0252969746878506, 0.023507498099856,
    0.881304771054552, 1.22455438983068, 0.879434379003448, 0.799071319536508,
    0.752627399778357, 0.68755580677248, 0.619393865502467, 1.03511573002578,
    0.402383261643063, 0.395061968298268, 0.3843251087703, 0.260852780265271,
    0.278178025770255, 0.359924428443045, 0.468267121413229, 0.410191603833112,
    0.241087622334221, 0.208878599816917, 0.386522274807038, 0.376516766580251,
    0.410435906289039, 0.441180670137384, 0.318441249000134, 0.388229279865325,
    0.241331924790148, 0.230350762806254, 0.428003898183351, 0.315024126750363,
    0.336010796961044, 0.274273854675231, 0.532686722514436, 0.336497845806299,
    0.394330616997087, 0.482419547124369, 0.36675556080939, 0.317952644088281,
    0.392377753416276, 0.276471020711968, 0.348943266459151, 0.371881244184046,
    0.344306187996348, 0.309657253052978, 0.347234705334266, 0.427516849338097,
    0.399210441849219, 0.528782551419412, 0.4238554246324, 0.428981108007056,
    0.2933076613051, 0.327958152315069, 0.407263086495194, 0.329178108528101,
    0.362851389714365, 0.407751691407047, 0.448746265938106, 0.475586858691426,
    0.496573528902107, 0.585882415242422, 0.465582906531237, 0.417755643567236,
    0.33210662586602, 0.394574919453014, 0.473391248721287, 0.461678735436213,
    0.374076854154185, 0.484616713161107, 0.531955371213256, 0.489985142925091,
    0.396281924511301, 0.430932415521269, 0.434105235315113, 0.366268511964136,
    0.621508560009497, 0.296478925032346, 0.544156489410183, 0.477295419816311,
    0.436789450197105, 0.493645011564189, 0.38920648968903, 0.482908152036222,
    0.311852863023117, 0.306240130803207, 0.310144301898232, 0.341379226725028,
    0.402138959187137, 0.379933888830021, 0.453870393246163, 0.469242775170336,
    0.46607151144309, 0.530491112544297, 0.60076774832134, 0.439229362623171,
    0.455090349459196, 0.458018866797114, 0.377736722793283, 0.357971564862234,
    0.424099727088326, 0.610773256548127, 0.336742148262225, 0.3843251087703,
    0.499989095085279, 0.65152352862326, 0.509994603312067, 0.369928380603234,
    0.459238823010147, 0.335277889593266, 0.448746265938106, 1.26180818025954,
    0.92091755844636, 1.06537500109547, 1.34599293929179, 1.11954634758056,
    0.95263953211841, 1.05805370775068, 1.06025087378741, 0.957519356970541,
    1.6200240475777, 0.722287213245558, 1.23179321164577, 1.30597401851784,
    1.12735468977061, 1.3074382771868, 1.22715768924957, 1.28132903573466,
    1.16078522456755, 1.12027769888174, 1.05317388289854, 0.92091755844636,
    1.03267581759972, 1.1446799352756, 0.959227918095426, 0.889195584774309,
    1.13174746577749, 1.21812939084648, 0.962887786734525, 1.28498890437376,
    1.22398486945572, 0.991438496679329, 1.16444509320665, 1.39796867580674,
    1.09392415497368, 1.01169070345563, 1.30036284236453, 1.2522912769446,
    1.5373019911478, 1.59562181118384, 1.17518195273461, 1.22154495702966,
    1.57463669703976, 1.22252061078676, 1.12223056246255, 1.47825081981058,
    1.71104149504289, 1.60391875828574, 1.55535858795396, 1.62587952618693,
    1.38430329894086, 1.51509692079068, 1.43847464542595, 1.69249629332488,
    1.34965280793089, 1.9574975471411, 2.71394820255137, 2.53313170776741,
    0.0130133849612163, 0.0131767719540332, -0.342518267474952,
    -0.564572083179304, -0.242390050076962, -0.300384652193969,
    -0.179758369497151, -0.184232060967137, -0.0619796886751415,
    -0.111433041234321, -0.152103953912942, -0.127701717519089,
    -0.161701772690986, -0.0781674494968004, -0.120218593248075,
    -0.184232060967137, -0.190658616017935, -0.215547901257043,
    0.00585547860923797, 0.00585547860923797, 0.00593795013894554,
    0.00593795013894554, 0.00553181675680068, 0.00601886560205487,
    -0.0676748924247591, -0.245072708892355, -0.14730504452392,
    -0.271346893403911, -0.315106598280071, -0.126318374313239,
    -0.384082362380972, -0.283954144982982, -0.114688332557873, -0.1241227643431,
    -0.122983723593176, -0.12485411564428, -0.231572275085885,
    -0.427516849338097, -0.645422191491498, -0.0788178853348715,
    -0.0600283811609283, -0.0558799076099773, -0.0552294717719062,
    -0.0565303434480484, -0.0570189483599008, -0.0582389045729337,
    -0.059133642866931, -0.0593779453228572, -0.0599474656978189,
    -0.0685696307187564, -0.0699514178580079, -0.0392860134061735,
    -0.0409945745310588, -0.13404580104018, -0.130468403930789,
    -0.0900417937080933, -0.0963049617660744, -0.109400818256998,
    -0.157145609691292, -0.153079607670049, -0.207820474530103,
    -0.147792093369175, -0.151451962008273, -0.10858699542611,
    -0.106064611470337, -0.0807691928490848, -0.0842672305619646,
    -0.0852428843190712, -0.0847542794072188, -0.140797574010013,
    -0.0902860961640195, -0.106147083000045, -0.091425136913943,
    -0.118673419116007, -0.19464370257607, -0.0815005441502652,
    -0.124447982262135, -0.0338366681790804, -0.0211485011368996,
    -0.031641058208941, -0.0283873229519873, -0.335441276586083,
    0.721800164400304, 0.514141520796419, 0.65933187081331, 0.663072654915517,
    0.577179334758375, 0.577829770596446, 0.431827153815266, -0.868208914563628,
    -0.833315677164332, -0.863329089711497, -0.99631832153146,
    -0.956299400757508, -0.880410032760555, -0.817941739173561,
    -0.687393975846261, -0.737172546324476, -0.728632852833246,
    -0.814037568078536, -0.820138905210299, -0.672752945223269,
    -0.742296673632534, -0.910912050219572, -0.828434296245602,
    -0.66494460303322, -0.991194194223403, -0.609308997879168,
    -0.728632852833246, -0.77914433067924, -0.867233260806522,
    -0.771578734878518, -0.669093076584171, -0.857960659947513,
    -0.98094438354069, -0.877481515422637, -0.875285905452497,
    -0.795736668816445, -0.829411506069307, -0.886754116281645,
    -0.768894519996527, -0.959227918095426, -0.975820256232632,
    -0.827947247400348, -0.920428953534507, -0.878458725246342,
    -0.779631379524494, -0.826482988731389, -0.766698910026387,
    -0.930677208150622, -0.798176581242511, -0.896759624508432,
    -0.673239994068523, -0.71423612466618, -0.854300791308415,
    -0.685196809809524, -0.995829716619608, -1.1097851418097, -0.946051146141393,
    -0.842099673111488, -1.00144244883952, -1.18640741717443, -0.803057962161241,
    -1.0858715103277, -0.708867694902196, -0.991925545524583, -0.81940599784252,
    -1.11246935669169, -0.943122628803475, -0.997782580200419,
    -0.978993076026476, -0.831607116039447, -0.81379326562261,
    -0.810377699439438, -0.824774427606504, -0.993145501737616,
    -0.910667747763646, -0.995585414163682, -0.701790704013327,
    -0.796225273728298, -0.97801742226937, -0.936534242826458,
    -0.997293975288567, -0.982408642209649, -0.943366931259401,
    -0.938974155252524, -0.884069901399653, -0.758401962924485,
    -0.888218374950604, -0.969476172711542, -0.763039041387289,
    -0.764503300056248, -0.803057962161241, -0.984848554635714,
    -0.845028190449406, -0.863084787255571, -0.85893631370462,
    -0.973380343806566, -0.883338550098473, -0.763039041387289,
    -0.893099755869334, -0.880898637672407, -0.934338632856319,
    -0.861864831042538, -0.92164890974754, -0.970451826468648,
    -0.953615185875516, -1.0858715103277, -2.45187391000647, -3.00237448293774,
    -2.91135703547254, -3.63754997587972, -2.62488206260039, -3.56873759877164,
    -1.94895785364987, -3.97160790925996, -3.92426769514121, -3.54189544995172,
    -3.05239891193848, -2.78690905321041, -2.44333421651524, -2.91038138171544,
    -3.49797235808274, -2.66490098337434, -2.75274716711229, -3.25786038950563,
    -3.3174001657547, -3.11803846925272, -2.75030569861963, -3.2685972490336,
    -2.81765537312535, -2.49091717702332, -2.78178492590235, -3.40988187188886,
    -3.713682090266, -2.77568358877059, -2.24031265137416, -2.59803991378047,
    -2.99822756545339, -3.81275218237718, -3.37937985442985, -2.72370940832223,
    -3.66878334463992, -3.26786434166582, -3.64779823049584, -3.53018138060005,
    -3.01042712758372, -3.90084266857106, -3.72319899358094, -4.18341347034819,
    -2.23689708519099, -4.11874956679106, -3.57166456004296, -2.56558658880724,
    -2.68369048754828, -3.0553258732098, -2.77836780365258, -2.33279435750832,
    -2.68027336529851, -3.31422734596086, -2.84596022454763, -2.75860264572153,
    -5.04893817002983, -3.96453091837109, -3.85813953291512, 0.697563649099165,
    0.688981163776483, 0.166180910394067, 4.08763445909233, 2.22087348939544,
    -2.68768024230621, -3.57499843272972, -1.50224303265579, -1.8411348869235,
    -3.43562232555722, -5.07724224341881, -4.68754792685123, -3.48821037427858,
    -2.49901883776714, -2.53833208233878, -4.27508446578466, -2.94108879996542,
    -1.46751084814941, 0.320225279355035, 0.318961753277251, 0.321244502976892,
    0.31949781822035, 0.319010769375096, 0.326206799358731, 0.376470862615602,
    -1.78648738405935, -1.61482056087313, -6.98846703812757, -9.23694826110927,
    -3.49382388453179, -5.19982061166314, -4.76717807895039, -4.86906387357111,
    -5.38157074640606, -5.6710971658774, -5.83747647476274, -6.30004605832646,
    -7.14967164754041, -6.39985139190531, -1.88933716193769, -2.36285133992023,
    -2.4663951235014, -2.43686409159969, -2.53781157806166, -2.574008021237,
    -2.56754645468774, -2.57717617283105, -2.64935976822425, -2.68024068789995,
    -2.76287716066694, -1.55337226894129, -1.5400694555928, -1.71865066070837,
    -2.36793656556334, -4.08137129103435, -4.2622002343511, -4.23095675115801,
    -4.2906614704665, -5.7470332158723, -6.50547174636184, -6.37572924749914,
    -5.74914635431273, -5.65406057072639, -3.20682296114942, -3.33101730455604,
    -3.36803612892856, -3.58350078062259, -3.74967780085017, -3.842803718556,
    -3.65500048474587, -3.78850866477635, -3.83215944499063, -4.02931930725608,
    -4.21505375052383, -4.9709737872241, -1.49646924754296, -0.440608037629226,
    -0.825593696670486, -0.508692953602625, -0.34617813611405,
    -0.393267823510478, 2.51153817158341, 40.9186236984312, 36.8092743299602,
    36.7823854991424, 34.2767930580309, 30.197329504646, 29.4625105090514,
    30.8039503974767, -10.6521581548247, -9.95964783200337, -10.2264743519393,
    -11.3993408775758, -10.5732204523618, -12.5151059833226, -10.4770765434886,
    -11.5807626822664, -10.9344854321792, -12.5904981880414, -9.96708505230973,
    -10.926195487377, -10.266734463036, -11.3146620673986, -10.7533063738778,
    -11.2823332277532, -11.4859694392006, -12.2560823593435, -10.2313549548248,
    -12.526445040624, -13.1214841297639, -11.8286923294332, -11.4650979179182,
    -9.79590839010196, -10.6829021406397, -9.30115779321971, -10.429619624375,
    -11.9774126165282, -10.5180773422861, -12.0029126579072, -10.336405010873,
    -11.9968105427421, -9.40292143661247, -10.3049265616236, -11.1926617778955,
    -10.2193981390838, -11.4635146201544, -10.3709294604886, -11.0274596333917,
    -11.7376632114685, -11.3820094078044, -11.110057982527, -10.5842086166453,
    -10.7189009633571, -9.76552930190421, -11.1365344556964, -10.5574816167537,
    -11.4744973382049, -10.3442164651963, -10.4878141810499, -11.3689594552781,
    -11.5158583664199, -11.034420697319, -10.9485234869954, -10.8947139259944,
    -9.97306657231343, -11.5982139691659, -10.9986607308245, -10.5814030285687,
    -11.5740568132613, -9.70269455463333, -10.7653822287136, -11.5075583071848,
    -10.6897418313724, -10.2967517657497, -12.1448057028023, -11.4897413446347,
    -10.2579465769224, -11.57234202787, -10.4652457691421, -10.5087961830608,
    -10.7589222182309, -11.1333600798359, -10.8529864440955, -11.9349592295613,
    -11.7671919092703, -11.1493439959332, -12.4375072761674, -11.1787388720076,
    -11.5585576119093, -11.8116043840844, -11.1704512613053, -11.4432600792782,
    -10.8909264598942, -11.550384372102, -11.6138275414128, -11.0398949396117,
    -11.4181272696165, -10.7215844002058, -11.4259332777067, -10.8736113288221,
    -11.5955266421507, -11.3565132565919, -11.5987017960444, -11.1509296277968,
    -10.8505434195363, -11.4370334787853, -12.5765908428195, -10.7757495224245,
    -4.00162132180712, -9.06629677137841, -9.8472733704771, -13.3626892351263,
    -12.3930924672534, -13.2541053518334, -10.2236749881291, -15.1373038310409,
    -18.0335320090119, -14.0171443932207, -15.411210453999, -13.4446814963216,
    -10.245021109724, -12.533527477746, -15.0892252633213, -12.0053580165663,
    -15.8093137544635, -13.9247871724144, -18.546689317685, -15.7126905770781,
    -14.0680168785175, -17.8804103875437, -13.4012462313312, -12.540238014951,
    -16.1240982469578, -16.5919989146871, -17.5867362725548, -14.1024300693712,
    -13.6928437754448, -15.6854858647408, -15.8287241292102, -15.3312939856457,
    -12.8999150288049, -15.7923604088755, -18.975556054797, -14.3832891978902,
    -17.6352917746867, -17.8302583610819, -16.7603785471845, -21.1277782199773,
    -22.2260943729246, -21.5496792249623, -20.0129879899543, -19.1894272942945,
    -19.2149320038732, -17.7377813231476, -15.9697644496095, -25.0059272132794,
    -18.5858562919964, -19.3779354262736, -23.4987522275405, -25.5764356905311,
    -20.4217145570852, -17.8518588995656, -24.5431568971245, -18.2269114074769,
    -7.15150391595986, 0.0, 0.0, 0.0, 4.0, 1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 13.0, 16.0, 17.0,
    14.0, 15.0, 16.0, 14.0, 7.0, 4.0, 8.0, 2.0, 5.0, 3.0, 5.0, 2.0, 3.0, 2.0,
    3.0, 5.0, 5.0, 4.0, 3.0, 3.0, 1.0, 2.0, 3.0, 1.0, 2.0, 3.0, 5.0, 3.0, 6.0,
    6.0, 3.0, 2.0, 5.0, 3.0, 7.0, 2.0, 6.0, 5.0, 3.0, 5.0, 5.0, 4.0, 4.0, 3.0,
    5.0, 5.0, 6.0, 4.0, 4.0, 4.0, 5.0, 5.0, 5.0, 7.0, 5.0, 6.0, 5.0, 4.0, 5.0,
    8.0, 3.0, 5.0, 6.0, 3.0, 8.0, 6.0, 5.0, 3.0, 4.0, 4.0, 6.0, 4.0, 3.0, 6.0,
    4.0, 3.0, 3.0, 5.0, 5.0, 4.0, 3.0, 3.0, 4.0, 3.0, 3.0, 7.0, 4.0, 6.0, 3.0,
    5.0, 4.0, 3.0, 7.0, 5.0, 7.0, 4.0, 3.0, 5.0, 7.0, 4.0, 6.0, 3.0, 4.0, 12.0,
    10.0, 11.0, 11.0, 9.0, 10.0, 8.0, 7.0, 11.0, 10.0, 8.0, 12.0, 10.0, 8.0, 9.0,
    10.0, 10.0, 9.0, 12.0, 8.0, 8.0, 8.0, 10.0, 9.0, 10.0, 6.0, 11.0, 10.0, 10.0,
    9.0, 10.0, 8.0, 10.0, 6.0, 9.0, 8.0, 11.0, 10.0, 8.0, 7.0, 9.0, 8.0, 6.0,
    9.0, 10.0, 11.0, 12.0, 8.0, 11.0, 10.0, 11.0, 7.0, 9.0, 9.0, 12.0, 11.0,
    16.0, 0.0, 0.0, 0.0, 18.0, 0.0, 11.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 43.0, 45.0, 47.0, 43.0, 46.0, 46.0,
    43.0, 36.0, 19.0, 38.0, 28.0, 26.0, 18.0, 39.0, 9.0, 38.0, 11.0, 41.0, 40.0,
    30.0, 39.0, 39.0, 39.0, 0.0, 24.0, 20.0, 0.0, 27.0, 11.0, 31.0, 29.0, 37.0,
    46.0, 38.0, 2.0, 29.0, 22.0, 40.0, 29.0, 30.0, 39.0, 19.0, 30.0, 21.0, 21.0,
    40.0, 39.0, 41.0, 30.0, 42.0, 20.0, 36.0, 29.0, 41.0, 27.0, 39.0, 39.0, 31.0,
    40.0, 35.0, 27.0, 37.0, 38.0, 19.0, 47.0, 38.0, 28.0, 39.0, 39.0, 30.0, 28.0,
    21.0, 27.0, 38.0, 37.0, 19.0, 40.0, 36.0, 37.0, 27.0, 40.0, 29.0, 37.0, 36.0,
    18.0, 36.0, 37.0, 37.0, 44.0, 20.0, 36.0, 37.0, 45.0, 27.0, 28.0, 28.0, 37.0,
    39.0, 38.0, 36.0, 29.0, 30.0, 28.0, 45.0, 27.0, 30.0, 43.0, 46.0, 44.0, 43.0,
    44.0, 41.0, 39.0, 40.0, 45.0, 41.0, 45.0, 41.0, 45.0, 42.0, 33.0, 45.0, 45.0,
    41.0, 41.0, 42.0, 38.0, 42.0, 41.0, 37.0, 42.0, 42.0, 46.0, 44.0, 40.0, 42.0,
    41.0, 40.0, 42.0, 37.0, 40.0, 36.0, 39.0, 42.0, 45.0, 36.0, 40.0, 43.0, 37.0,
    38.0, 44.0, 43.0, 42.0, 42.0, 40.0, 41.0, 41.0, 41.0, 47.0, 41.0, 39.0, 42.0,
    46.0, 0.0, 0.0, 0.0, 4.5, 0.0, 5.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.5, 3.30769230769231, 2.8125,
    2.76470588235294, 3.07142857142857, 3.06666666666667, 2.875,
    3.07142857142857, 5.14285714285714, 4.75, 4.75, 14.0, 5.2, 6.0, 7.8, 4.5,
    12.6666666666667, 5.5, 13.6666666666667, 8.0, 6.0, 9.75, 13.0, 13.0, 0.0,
    12.0, 6.66666666666667, 0.0, 13.5, 3.66666666666667, 6.2, 9.66666666666667,
    6.16666666666667, 7.66666666666667, 12.6666666666667, 1.0, 5.8,
    7.33333333333333, 5.71428571428571, 14.5, 5.0, 7.8, 6.33333333333333, 6.0,
    4.2, 5.25, 10.0, 13.0, 8.2, 6.0, 7.0, 5.0, 9.0, 7.25, 8.2, 5.4, 7.8,
    5.57142857142857, 6.2, 6.66666666666667, 7.0, 6.75, 7.4, 4.75,
    6.33333333333333, 9.4, 6.33333333333333, 9.33333333333333, 4.875, 6.5, 6.0,
    9.33333333333333, 5.25, 6.75, 6.33333333333333, 9.25, 6.33333333333333,
    6.66666666666667, 9.0, 12.3333333333333, 9.0, 8.0, 5.8, 9.25, 12.0, 6.0, 9.0,
    12.3333333333333, 12.3333333333333, 6.28571428571429, 5.0, 6.0,
    12.3333333333333, 9.0, 6.75, 9.33333333333333, 4.0, 7.4, 5.57142857142857,
    9.5, 12.0, 5.8, 4.28571428571429, 7.0, 7.5, 9.0, 7.5, 3.58333333333333, 4.6,
    4.0, 3.90909090909091, 4.88888888888889, 4.1, 4.875, 5.71428571428571,
    4.09090909090909, 4.1, 5.625, 3.41666666666667, 4.5, 5.25, 3.66666666666667,
    4.5, 4.5, 4.55555555555556, 3.41666666666667, 5.25, 4.75, 5.25, 4.1,
    4.11111111111111, 4.2, 7.0, 4.18181818181818, 4.4, 4.0, 4.66666666666667,
    4.1, 5.0, 4.2, 6.16666666666667, 4.44444444444445, 4.5, 3.54545454545455,
    4.2, 5.625, 5.14285714285714, 4.44444444444445, 5.375, 6.16666666666667,
    4.22222222222222, 4.4, 3.90909090909091, 3.5, 5.25, 3.63636363636364, 4.1,
    3.72727272727273, 5.85714285714286, 5.22222222222222, 4.55555555555556, 3.25,
    3.81818181818182, 2.875, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 2.0, 0.0, 1.0, 1.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 1.0, 1.0, 0.0, 0.0,
    1.0, 0.0, 0.0, 1.0, 1.0, 2.0, 4.0, 0.0, 1.0, 2.0, 1.0, 0.0, 3.0, 4.0, 1.0,
    3.0, 3.0, 2.0, 1.0, 3.0, 4.0, 2.0, 3.0, 5.0, 9.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 6.0, 49.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 2.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 2.0, 1.0, 4.0, 2.0,
    2.0, 2.0, 0.0, 0.0, 3.0, 1.0, 0.0, 1.0, 2.0, 2.0, 0.0, 1.0, 3.0, 0.0, 2.0,
    0.0, 1.0, 1.0, 1.0, 2.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0,
    1.0, 3.0, 2.0, 1.0, 1.0, 2.0, 0.0, 0.0, 2.0, 1.0, 0.0, 0.0, 2.0, 2.0, 2.0,
    0.0, 1.0, 0.0, 1.0, 10.0, 11.0, 9.0, 13.0, 9.0, 9.0, 7.0, 7.0, 11.0, 10.0,
    7.0, 8.0, 6.0, 8.0, 6.0, 7.0, 7.0, 10.0, 13.0, 8.0, 8.0, 7.0, 14.0, 9.0, 6.0,
    7.0, 10.0, 11.0, 7.0, 8.0, 10.0, 12.0, 13.0, 6.0, 11.0, 9.0, 10.0, 10.0, 2.0,
    10.0, 6.0, 7.0, 7.0, 7.0, 6.0, 8.0, 11.0, 6.0, 10.0, 6.0, 4.0, 5.0, 4.0, 6.0,
    8.0, 6.0, 10.0, 50.0, 50.0, 50.0, 49.0, 50.0, 47.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 44.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    49.0, 50.0, 50.0, 50.0, 50.0, 50.0, 48.0, 49.0, 50.0, 50.0, 50.0, 50.0, 49.0,
    49.0, 50.0, 50.0, 50.0, 50.0, 49.0, 50.0, 50.0, 50.0, 48.0, 50.0, 50.0, 49.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 49.0, 50.0, 49.0, 49.0,
    49.0, 50.0, 50.0, 49.0, 50.0, 50.0, 49.0, 48.0, 49.0, 46.0, 48.0, 48.0, 48.0,
    50.0, 50.0, 47.0, 49.0, 50.0, 49.0, 48.0, 48.0, 50.0, 49.0, 47.0, 50.0, 48.0,
    50.0, 49.0, 49.0, 49.0, 48.0, 50.0, 49.0, 50.0, 50.0, 50.0, 50.0, 49.0, 50.0,
    49.0, 49.0, 49.0, 47.0, 48.0, 49.0, 49.0, 48.0, 50.0, 50.0, 48.0, 49.0, 50.0,
    50.0, 48.0, 48.0, 48.0, 50.0, 49.0, 50.0, 49.0, 39.0, 39.0, 41.0, 36.0, 41.0,
    41.0, 43.0, 43.0, 39.0, 39.0, 43.0, 41.0, 42.0, 42.0, 43.0, 42.0, 42.0, 40.0,
    37.0, 42.0, 42.0, 43.0, 36.0, 41.0, 44.0, 43.0, 38.0, 39.0, 42.0, 41.0, 40.0,
    38.0, 36.0, 44.0, 39.0, 40.0, 39.0, 38.0, 44.0, 40.0, 43.0, 41.0, 42.0, 43.0,
    41.0, 38.0, 38.0, 41.0, 37.0, 42.0, 45.0, 42.0, 42.0, 42.0, 39.0, 39.0, 31.0,
    0.000325217919035541, 0.000244302455926215, 0.0200094603869761,
    0.00618225259487176, 0.0644196011012073, 0.0190322505632713,
    0.0359513626861106, 0.0663724646820187, 0.0123629491231453,
    0.0156166843800989, 0.00171011719148351, 0.00634563958768863,
    0.0088664674768637, 0.00252082788917501, 0.0216371060487521,
    0.00593795013894555, 0.0361956651420368, 0.0846749200107078,
    0.0001633869928169, 0.000407689448743118, 0.000407689448743118,
    1.55606659825447E-6, 7.93593965110681E-5, 0.000488604911852441,
    0.0069136038960522, 0.0353818423111488, 0.019359024548905,
    0.0457934839200809, 0.0211485011368996, 0.00431186054376787,
    0.0165938942038038, 0.00772742672694027, 0.00463707846280337,
    0.000894738293997333, 0.000406133382144791, 0.000406133382144902,
    0.0579961581836057, 0.00439277600687715, 0.120056762321856,
    0.0180581528727629, 0.00203377904392077, 0.000650435838071006,
    0.000487048845254301, 0.00049016097845056, 0.000488604911852431,
    0.00203377904392077, 0.00097565375710662, 0.00048704884525419,
    0.00113904074992344, 0.00406755808784165, 0.00528751430087437,
    0.00276668641169942, 0.00138334320584976, 0.00203377904392077,
    0.0097596497042628, 0.000323661852437374, 0.00536687369738553,
    0.0307463199149436, 0.00910921386619157, 0.0627935115060295,
    0.0473402141187475, 0.00146581473555729, 0.000812266764289693,
    0.00252394002237144, 0.00772587066034186, 0.00447213540338826,
    0.00227963756644528, 0.00211469450703011, 0.00292696127131986,
    0.00211469450703011, 0.00130398380933866, 0.00276513034510123,
    0.0029285173379181, 0.0152914664610634, 0.0253778901509599,
    0.0179772374096536, 0.108017475051149, 0.00284760187480881,
    0.00138334320584965, 0.00821447557219435, 0.00756403973412323,
    0.0712522895341501, 0.0710889025413332, 0.145839229788363,
    0.0220432394308969, 0.0100879797564948, 0.00984367730056853,
    0.0414831794429112, 0.030013412547165, 0.0929703110460118,
    0.0692994259533386, 0.0897974912521669, 0.176668021233014,
    0.0195208554751236, 0.107121180690553, 0.406043130282161, 0.0780849779670928,
    0.00902829840308228, 0.163002644367129, 0.169835332800071, 0.200580096648417,
    0.27842388429278, 0.169348283954817, 0.0222066264237139, 0.166906815462153,
    0.135916193091283, 0.0568571174336822, 0.22425098174109, 0.0175711040275085,
    0.0275735001210994, 0.0561226539993052, 0.0788178853348716,
    0.183500709665957, 0.243284788370959, 0.277445118402477, 0.0253778901509598,
    0.130792065783226, 0.283302153078313, 0.0324548810398289, 0.229376665115745,
    0.561481734915167, 0.147873008832284, 0.306484433259133, 0.124692284718062,
    0.264756951360296, 0.0631996448881743, 0.259390077662911, 0.14445899871571,
    0.0358704472230011, 0.0224493728130418, 0.0680810258069041,
    0.00414847355095094, 0.131769275606931, 0.059784078705002, 0.14738595998703,
    0.405065920458457, 0.217173990852221, 0.00732129334479525, 0.206192828868327,
    0.055635605154051, 0.153972789897448, 0.192284705613113, 0.0173252455049844,
    0.0817448466061914, 0.23620779748209, 0.134451934422324, 0.199848745347237,
    0.0102482546161153, 0.0636882498000269, 0.357238657494455, 0.120789669689635,
    0.18472066587899, 0.00732129334479525, 0.33503358713734, 0.21351567827972,
    0.00975964970426269, 0.142017530223046, 0.0268421488199191,
    0.215221127271409, 0.0712538456007483, 0.122008069836069, 0.0358704472230011,
    0.415802779986424, 0.0932146135019375, 0.0383103596490669,
    0.00610133713176242, 0.0373347058919604, 0.220102508190139,
    0.312584214324298, 0.366268511964136, 0.12469384078466, 0.0300149686137634,
    0.0485586142651819, 0.291600656246813, 0.0114697668957462,
    0.0275735001210995, 0.0556371612206494, 0.223762376829237, 0.224006679285163,
    0.14592170131807, 0.046851609206895, 0.107611341669004, 0.28696357778401,
    0.134453490488923, 0.285742065504379, 0.130549319393898, 0.180572192328039,
    0.00414847355095094, 0.895783970751326, 1.41675662391409, 1.4858148595447,
    1.21861643969174, 0.827702944944422, 1.59928167982294, 1.20714822886259,
    1.27107922505195, 0.762063387630182, 1.29938563254082, 1.28108473327873,
    1.47434664871555, 1.53168925892789, 0.927015783444925, 1.18177033871163,
    1.03048020762958, 1.1971442767024, 1.22715768924957, 0.794759458992741,
    1.19787562800358, 1.06903486973457, 1.37600790790555, 1.09587701855449,
    1.52363817034851, 0.936045637914606, 0.77987568198042, 1.65076881142604,
    1.41480531639988, 1.03292012005564, 1.05634514662579, 1.29084593904959,
    1.18323459738059, 1.56341123259994, 1.07489034834381, 1.14345997906256,
    1.33184051358065, 1.20665962395074, 0.897003926964358, 0.467291467656122,
    0.860400572373578, 0.318439692933536, 1.27181213241972, 0.733512677685378,
    1.0068108786035, 0.651279226167334, 0.596374972314463, 0.519021345648551,
    0.0824777539739703, 1.18030608004267, 1.08709146654073, 0.81867464654134,
    0.910180698918392, 0.156414258390112, 1.05048966801655, 1.53046930271486,
    0.773287296003403, 0.326249591190183, -1.91304658256201E-5,
    7.18536635077104E-6, 0.0200094603869761, 0.00206075086495725,
    -0.0644196011012073, 0.00271889293761018, 0.00719027253722212,
    -0.0663724646820187, 0.00618147456157267, -0.0156166843800989,
    0.000190013021275945, 0.00019830123711527, 0.00221661686921593,
    -0.00252082788917501, -0.0216371060487521, -0.00593795013894555,
    0.0361956651420368, -0.0211687300026769, 2.04233741021125E-5,
    1.69870603642966E-5, -0.000407689448743118, 2.22295228322067E-7,
    5.29062643407121E-6, 1.74501754233015E-5, 0.00138272077921044,
    -0.0353818423111488, -0.0019359024548905, 0.00269373434824005,
    0.0105742505684498, -0.000159698538658069, 0.00207423677547547,
    -0.00193185668173507, 0.000144908701962605, -0.000127819756285333,
    -1.84606082793087E-5, 0.000406133382144902, -0.0289980790918029,
    0.000549097000859644, 0.120056762321856, -0.00902907643638146,
    0.000169481586993397, -5.91305307337278E-5, -4.87048845254301E-5,
    0.00012254024461264, 0.000488604911852431, 0.00203377904392077,
    -0.00097565375710662, 0.00048704884525419, -4.55616299969375E-5,
    0.000193693242278174, 0.00176250476695812, -0.000307409601299935,
    -0.00138334320584976, -0.000169481586993397, -0.0048798248521314,
    0.000107887284145791, 0.00536687369738553, 0.0153731599574718,
    -0.00910921386619157, -0.0627935115060295, 0.0236701070593737,
    0.000732907367778646, 0.000812266764289693, -0.000157746251398215,
    0.000858430073371318, 0.000406557763944387, 0.00056990939161132,
    0.000105734725351506, 0.000225150867024605, -0.00211469450703011,
    0.00065199190466933, -0.00276513034510123, 0.0029285173379181,
    -0.0152914664610634, -0.00507557803019197, 0.0179772374096536,
    -0.0180029125085248, -0.000113904074992353, 0.00138334320584965,
    -0.000186692626640781, 0.000756403973412323, 0.035626144767075,
    -0.00444305640883332, -0.0364598074470908, 0.00183693661924141,
    0.00201759595129895, 0.000307614915642767, 0.00218332523383743,
    -0.00157965329195605, 0.010330034560668, 0.0692994259533386,
    0.0897974912521669, 0.0135898477871549, 0.0195208554751236,
    0.00357070602301843, -0.0812086260564323, -0.00205486784123928,
    0.00902829840308228, 0.00776203068414898, 0.169835332800071,
    0.200580096648417, 0.027842388429278, -0.0169348283954817,
    0.0222066264237139, 0.00575540742972943, -0.00503393307745494,
    -0.00138675896179713, 0.0203864528855536, -0.000976172445972694,
    0.0275735001210994, 0.0561226539993052, -0.0157635770669743,
    0.00873812903171225, -0.0486569576741919, -0.0115602132667699,
    0.0253778901509598, -0.00769365092842505, 0.283302153078313,
    0.00147522186544677, 0.229376665115745, -0.0207956198116729,
    -0.00821516715734913, 0.0145944968218635, -0.00337006174913681,
    0.264756951360296, -0.0126399289776349, -0.0103756031065164,
    0.14445899871571, 0.0358704472230011, 0.0224493728130418,
    -0.00219616212280336, -0.000218340713207944, 0.131769275606931,
    0.059784078705002, 0.14738595998703, 0.405065920458457, 0.217173990852221,
    -0.000430664314399721, 0.206192828868327, 0.00264931453114528,
    -0.0219961128424926, 0.192284705613113, -0.000481256819582899,
    0.0817448466061914, 0.23620779748209, 0.134451934422324, -0.0111027080748465,
    0.0102482546161153, -0.00172130404864937, 0.357238657494455,
    -0.00525172476911457, 0.18472066587899, -0.000406738519155292,
    -0.033503358713734, 0.0152511198771229, 0.000487982485213134,
    -0.0355043825557614, 0.00134210744099595, -0.00827773566428497,
    0.00356269228003742, 0.122008069836069, 0.0358704472230011,
    0.415802779986424, 0.00443879111913988, 0.0383103596490669,
    0.000305066856588121, 0.0019649845206295, 0.220102508190139,
    0.0104194738108099, 0.366268511964136, -0.00733493181086234,
    -0.00176558638904491, 0.0485586142651819, 0.0243000546872344,
    0.0114697668957462, -0.00459558335351658, -0.00120950350479673,
    0.223762376829237, 0.224006679285163, 0.14592170131807, -0.00520573435632167,
    -0.0179352236115006, 0.28696357778401, 0.134453490488923, 0.285742065504379,
    0.00652746596969489, 0.180572192328039, 0.00414847355095094,
    -0.052693174750078, 1.41675662391409, 0.106129632824621, 1.21861643969174,
    0.827702944944422, 1.59928167982294, 1.20714822886259, 1.27107922505195,
    -0.0423368548683434, 1.29938563254082, 1.28108473327873, 1.47434664871555,
    1.53168925892789, 0.066215413103209, 1.18177033871163, 1.03048020762958,
    1.1971442767024, 1.22715768924957, 0.794759458992741, 1.19787562800358,
    1.06903486973457, 1.37600790790555, -0.0438350807421795, 1.52363817034851,
    0.936045637914606, 0.77987568198042, 1.65076881142604, 1.41480531639988,
    1.03292012005564, 1.05634514662579, 1.29084593904959, 1.18323459738059,
    -0.0422543576378362, 1.07489034834381, 1.14345997906256, 1.33184051358065,
    -0.603329811975369, 0.897003926964358, 0.467291467656122, 0.860400572373578,
    0.0122476804974437, 0.0489158512469125, 0.733512677685378, 1.0068108786035,
    0.651279226167334, 0.596374972314463, 0.519021345648551, 0.0824777539739703,
    1.18030608004267, 1.08709146654073, 0.81867464654134, -0.0267700205564233,
    0.156414258390112, -0.0456734638268066, 1.53046930271486, -0.386643648001702,
    0.326249591190183, 1.09830915064756, 1.09660525772247, 1.10856674166327,
    20.4780745112407, 29.2302940223938, 25.7006806061021, 23.0668527166284,
    22.7573541823684, 20.482806509766, 20.3767341179627, 20.5663330526272,
    19.4767207581973, 19.7984157424544, 19.4418290768646, 19.1234656311944,
    18.8645657145099, 20.9478090034561, 9.24878370365561, 1.0278753521441,
    1.02730894390233, 1.02763571788796, 1.02372843465974, 1.0187676943445,
    1.03568525040074, 1.94977323254736, 32.7799168990725, 33.6618767741649,
    31.0102241418081, 29.5806175199256, 31.0783938634111, 32.0486760786206,
    30.5709200879882, 31.5906618800232, 31.159888183295, 31.1000278573267,
    31.1027945437384, 30.4111711788781, 40.0459869963622, 35.6732935850022,
    33.0522363341003, 32.9730963429796, 32.5075865715785, 32.5550139254267,
    32.6322041650998, 32.6623731843068, 32.6599410522137, 32.2392568913086,
    32.0626744537385, 31.9982501844375, 31.7991281221916, 31.928703343895,
    31.5586209126985, 30.9706969380792, 30.4421618012489, 27.0365933459704,
    27.1422409315924, 27.1214285408408, 27.1533092333058, 28.3123003130856,
    28.9895316179786, 29.4568215295682, 27.9681559560154, 27.5016860915232,
    26.4180086350991, 24.7813082217878, 24.1109158295941, 23.9834553023978,
    23.8835760556555, 23.6185701336395, 24.0445511452451, 24.2212923016082,
    23.7658736220968, 23.650940987083, 23.5630963594117, 23.6718405175642,
    23.9749996365028, 29.9372928813115, 31.1189791924268, 31.0927050079153,
    34.1741416786438, 35.9809995305409, 25.1925595069744, 0.416291384898277,
    1.49566009291187, 1.8158861503002, 2.51109002440311, 2.22810686497744,
    2.45780874801222, 2.33214547773685, -49.4362124855995, -48.5091982582212,
    -48.9154934712923, -48.2917706285135, -48.7866464887569, -47.5385021292966,
    -48.5475117300035, -48.5797223085874, -48.2324860471865, -48.1861121504253,
    -48.2046698006761, -48.0272673160087, -48.718557682617, -48.5284779233736,
    -48.0170143931928, -48.6107067066919, -46.9033281924213, -47.7446980702982,
    -48.6024144277897, -46.1197957539349, -46.663213111711, -48.5160231663212,
    -48.5392101147018, -48.4523426968541, -48.4208696938378, -49.1775537592378,
    -47.8862269956094, -48.903770065541, -48.3469269651553, -48.7378357917028,
    -47.669541609669, -49.3095704494339, -49.3122484400495, -47.8269268536164,
    -47.9640661150539, -48.8698618182984, -48.6331560795049, -48.39549491582,
    -48.1068134405116, -48.0484920644089, -48.4538069555231, -48.32862295376,
    -48.5101754680449, -48.4425830471499, -48.580452103822, -48.5252973232468,
    -48.4352586416719, -48.3842570028474, -49.9735362865768, -48.2588333668282,
    -48.4032939216105, -48.3283770952374, -49.8464103136991, -47.8491319239735,
    -47.9713905205319, -48.696600026849, -48.6878144748352, -48.5941097003549,
    -47.5411863441786, -47.8369323618432, -49.8085870028952, -47.6500145299275,
    -47.6402533241566, -48.8996215919901, -49.3334825248493, -47.6858834210839,
    -48.323746241041, -49.0235871936067, -48.0350803263985, -48.4272059970259,
    -49.1777996177603, -47.4940951007156, -49.2675971090125, -48.5357898803188,
    -47.8274185706614, -48.2939740188166, -48.8947433232045, -48.4245280064103,
    -47.4953135008621, -48.4840631144596, -47.8889080983582, -48.8596011151495,
    -47.7603163107449, -49.7836946055229, -46.9958098985554, -49.2812609298118,
    -47.4665169323947, -49.1953722778544, -47.7615362669579, -49.043347683338,
    -47.3618341080636, -49.0282196038698, -47.1924920483752, -48.6717138537431,
    -47.4799380068047, -49.3098147518898, -47.9059921535404, -47.9191642572946,
    -49.518936098096, -44.7181609852227, -54.6596292911586, -49.6553393400326,
    -51.3305119439191, -50.4693737959779, -47.3579314930352, -50.0174593784458,
    -53.5420358071589, -50.750732421875, -49.5152762294569, -50.9720533302116,
    -53.5098190043086, -48.7327256689941, -47.7068732034277, -53.4190489714325,
    -48.333263144356, -52.580369532704, -42.8682747797506, -44.4546286582752,
    -53.9734521593922, -46.7015250274266, -53.346817916008, -46.794501562739,
    -47.0580354457531, -52.4110228048158, -47.9250228480371, -54.3402077202016,
    -46.3779596511184, -49.309817864023, -43.9512162166738, -50.815877150011,
    -51.9956992933416, -49.1526722543317, -51.2412061697119, -49.7849176738691,
    -50.5948067683968, -40.8863749962156, -50.0740721934235, -43.9600017686876,
    -41.4388315448609, -47.8054609148935, -50.3988637502112, -47.2683767481724,
    -48.7756622146398, -51.474729752312, -49.6526629054836, -43.4700166257625,
    -50.2565972493324, -42.8633887306321, -48.2053995959107, -48.4530740481553,
    -51.2792768951049, -47.771292804529, -54.4263468988813, -51.135306501301,
    -45.2647589431256, -51.1980159852108, 0.000493271482748432,
    0.000429727678110511, 0.0100788159098499, 0.159063283251629,
    0.0753027085990925, 0.0866877031902511, 0.0443883282658408,
    0.0268893678462579, 0.0297709266314517, 0.0142955007143031,
    0.0129318880932793, 0.0108899920957711, 0.0223286627577795,
    0.0172896953884178, 0.0168989819588085, 0.0111579528267867,
    0.040904377151838, 0.170078541237914, 0.000630715223574303,
    0.000531517405552706, 0.000532814131096944, 0.000513484005839969,
    0.000479425001732981, 0.000544648631403144, 0.0737688084092456,
    0.0636076111133392, 0.0159794331550905, 0.0433971765459336,
    0.0240329576802638, 0.00899065706927205, 0.0508847781452208,
    0.0295853374801121, 0.00803838353690654, 0.00252726652527175,
    0.00236544636964573, 0.0016520109396459, 0.0174992813579882,
    0.203753090906433, 0.164840928456815, 0.0157626154462057,
    0.00628189361333325, 0.0022476701356885, 0.00230367212509002,
    0.00128203224645846, 0.00153771096583796, 0.00191599108214474,
    0.00252678271970046, 0.00151605352234738, 0.00175215524283125,
    0.00270242939694979, 0.00560569984925872, 0.00597568880287484,
    0.00310103383238878, 0.0283737896648701, 0.0160135114810882,
    0.00182657922001304, 0.00293733410269551, 0.00646054288393522,
    0.01806028172263, 0.0112195899163295, 0.0349295374671393, 0.0135008955675335,
    0.00961150408005823, 0.0253949432093433, 0.0118128784703045,
    0.00596886991365269, 0.00537952869060773, 0.00226454529255342,
    0.00279983731924674, 0.013251600886599, 0.00436708159819045,
    0.00372947010234487, 0.0027041990147916, 0.0043259943374733,
    0.0184763332409848, 0.0166032332509355, 0.0489342245981618,
    0.00384252814357936, 0.00431009755547073, 0.0289078300275955,
    0.0148160324157793, 0.154685669484906, 0.0928466660014032,
    0.0905422539993967, 0.0215985660570258, 0.0204538803058227,
    0.0215832077149024, 0.0202299143203869, 0.0731935461699001,
    0.374355712467968, 0.303976564197611, 0.354869575119358, 0.305375237862015,
    0.329790796416842, 0.269692950250347, 0.248406361930464, 0.262447559392895,
    0.227274507543658, 0.216942323465736, 0.201068643436984, 0.2366108015694,
    0.247626257755785, 0.248341462565028, 0.234604958711659, 0.278193663404207,
    0.29884146103437, 0.231777190280674, 0.241164849231949, 0.28976946392638,
    0.293686053502677, 0.247202649577991, 0.263937530784086, 0.261929218669259,
    0.296800491744714, 0.305112815942661, 0.306972596404449, 0.286492420765054,
    0.310831388681824, 0.340720775108207, 0.301437623628719, 0.365573964444872,
    0.400801046212304, 0.327885641516413, 0.366865765119219, 0.305620703039304,
    0.28301944095654, 0.271062325869223, 0.317259893728255, 0.266387828110991,
    0.274377017256057, 0.241107852765802, 0.249410772917007, 0.221856889086626,
    0.260156913824413, 0.303710900852123, 0.300851906523479, 0.358667427524742,
    0.371734941870464, 0.313077770801349, 0.272206949853093, 0.335980890924732,
    0.371520562342764, 0.31385168296282, 0.330779292702497, 0.400588058782691,
    0.328464699688889, 0.284746001262126, 0.322210426820858, 0.285004127469981,
    0.344287425682345, 0.309002356320445, 0.319018442282469, 0.286893375559221,
    0.29909563790085, 0.265334425269188, 0.334591317235181, 0.35171931349734,
    0.334505188349297, 0.265933246693477, 0.312205293017022, 0.359745832053921,
    0.354353077600241, 0.314235895744068, 0.298379525551437, 0.286105467935764,
    0.335295352147034, 0.359992621712014, 0.365397504587379, 0.282688461900131,
    0.315648394488173, 0.36804673696619, 0.314227900404483, 0.369330544669486,
    0.31079957331706, 0.305137577991432, 0.330016951684355, 0.331769314037209,
    0.405148690622708, 0.346847443382788, 0.372422482545584, 0.329458271395646,
    0.346484938989179, 0.358747881892096, 0.370025662771518, 0.34862682563292,
    0.328439378032277, 0.331832751459967, 0.352237862844777, 1.35227596007995,
    1.41107580171184, 1.52695060104623, 1.65514220305154, 1.51498174791729,
    1.43199029900463, 1.62694800114228, 1.67902412144008, 1.62904106620513,
    1.7199819802878, 1.65375341842661, 1.44175071034126, 1.4034603649186,
    1.43532607479116, 1.4439773158811, 1.49776402148953, 1.64076651155483,
    1.30131806228438, 1.25276207269118, 1.58749815726035, 1.54253069543398,
    1.74874510218114, 1.52648926302464, 1.47744139818509, 1.71501960063327,
    1.49249178425003, 1.70116719930083, 1.40243565916934, 1.57086922123467,
    1.46408883938774, 1.61589561500364, 1.6507448428491, 1.5868267629766,
    1.55131839466935, 1.51657465623174, 1.65477808179108, 1.3599475357256,
    1.68143804750907, 1.52263745008754, 1.4652495656898, 1.52479215133363,
    1.61115448213242, 1.52341567568286, 1.62540544650323, 1.60307457500748,
    1.70757875931751, 1.37333765072839, 1.80876580169346, 1.44385480973672,
    1.59104066996648, 1.7310190485754, 1.6796501645198, 1.8005707311808,
    1.92247847057831, 1.75505342675507, 1.37470226391613, 1.92701817836325,
    0.0233441111070391, 0.0229379777248942, 0.0561242100659035,
    0.680153597964575, 0.884477590848397, 0.681292638714499, 0.547085006748101,
    0.572299509906244, 0.458832689628002, 0.429306325926091, 0.43508088907222,
    0.409377781002224, 0.475099809846172, 0.421009378824189, 0.436545147741179,
    0.4174304256482, 0.542773146204333, 0.424018811625217, 0.0222875418868231,
    0.02187985243808, 0.0215546345190445, 0.0216355499821538, 0.0213912475262276,
    0.02187985243808, 0.432233287197411, 0.797689532397257, 0.698619440286081,
    0.729364204134427, 0.65014018541741, 0.639241494963224, 0.747420800940591,
    0.668115866760466, 0.656241522549173, 0.627609897141259, 0.627609897141259,
    0.625006597722377, 0.622404854370092, 1.14280954322449, 0.979154906952695,
    0.698698799682592, 0.674459950281556, 0.653880969519618, 0.65713626084317,
    0.654696348417104, 0.656159051019465, 0.658030999137167, 0.652010577468514,
    0.644446537734391, 0.643715186433211, 0.644609924727208, 0.648595011285342,
    0.641518020396473, 0.625983807546082, 0.684221156052417, 0.568964859186181,
    0.547327753137429, 0.553916139114445, 0.560748827547388, 0.596457443844171,
    0.599384405115491, 0.658762350438348, 0.607681352217393, 0.56904577464929,
    0.588649101654121, 0.526750328442089, 0.493806842490408, 0.49218075289523,
    0.483558587874293, 0.478271073573418, 0.513735387414275, 0.493481624571372,
    0.48518467746947, 0.479816247705487, 0.492505970814266, 0.548874483336095,
    0.506821783518222, 0.680805589869245, 0.629724591648289, 0.630862076331615,
    0.743110496463422, 0.745306106433561, 0.716676037092246, 0.181792148541072,
    0.263536995147263, 0.112898855969878, 0.0966301796851099, 0.0929703110460114,
    0.0865437559952134, 0.304369738752103, -0.573925599501421,
    -0.536103844764207, -0.597352182138168, -0.516825735678411,
    -0.471926990052328, -0.530491112544297, -0.431421020433121,
    -0.536103844764207, -0.586126717698348, -0.61492017403248, -0.5866153226102,
    -0.562945993584126, -0.597107879682242, -0.6261456384723, -0.533419629882215,
    -0.549769221630093, -0.430688113065343, -0.543912186954256,
    -0.507553134819403, -0.507553134819403, -0.515850081921305,
    -0.463142994105172, -0.578562677964225, -0.53732380097724,
    -0.512434515738132, -0.541227972072265, -0.56343304242938,
    -0.508041739731255, -0.544156489410183, -0.42946815685231,
    -0.450453270996393, -0.499257743784099, -0.533906678727469,
    -0.507310388430075, -0.482908152036222, -0.493645011564189,
    -0.495840621534329, -0.496084923990255, -0.522682770354248,
    -0.548547709350462, -0.477295419816311, -0.554160441570372,
    -0.625169984715194, -0.594910713645504, -0.569534379561143,
    -0.584662459029389, -0.595887923469209, -0.532686722514436,
    -0.536348147220133, -0.581733941691471, -0.56904577464929,
    -0.449477617239286, -0.433128025491408, -0.536348147220133,
    -0.509505998400214, -0.447769056114401, -0.481199590911336,
    -0.561968783760421, -0.553429090269191, -0.528051200118231,
    -0.613700217819447, -0.498280533960394, -0.519998555472256,
    -0.612724564062341, -0.579049726809479, -0.589543839948119,
    -0.324785332521224, -0.399699046761071, -0.556113305151183,
    -0.568557169737438, -0.56343304242938, -0.547572055593355,
    -0.539763713403305, -0.542447928285297, -0.540252318315158,
    -0.483639503337402, -0.453137485878384, -0.470707033839295,
    -0.45265043703313, -0.536590893609461, -0.535370937396428,
    -0.497060577747361, -0.517314340590264, -0.48095528845541,
    -0.479003980941197, -0.489252235557312, -0.521462814141215,
    -0.528539805030084, -0.42092846336108, -0.453137485878384,
    -0.510237349701395, -0.513410169495239, -0.498769138872247,
    -0.479979634698303, -0.519998555472256, -0.524391331479133,
    -0.540496620771084, -0.506088876150444, -0.432396674190228, 1.38625616252167,
    1.18835872468865, 1.33184051358065, 1.82109274913796, 1.08245594414453,
    1.28620886058679, 1.26595665381049, 1.37332369302356, 1.31915079047187,
    1.30621832097377, 1.41212110151788, 1.25033841336379, 0.868453217019555,
    1.46263257936388, 1.34013746068255, 1.44140316276387, 1.29133298789485,
    1.28767311925575, 1.37234648319986, 1.0785517730495, 1.29694572011476,
    1.40480136423969, 1.52973795141368, 1.23008620658749, 1.48361769350796,
    1.40919414024656, 1.35648549636383, 1.33354907470554, 1.34477298307876,
    1.4306663032359, 1.3635624872527, 1.3672223558918, 1.43481477678685,
    1.31915079047187, 1.32964490361051, 1.41285400888566, 1.72226695948272,
    1.48532625463285, 1.40797418403353, 1.46458544294469, 1.51680392584897,
    1.27278778617683, 1.29401720277684, 1.48093503469257, 1.81840853425597,
    1.44042750900676, 1.7442277273839, 1.49240324552172, 1.61343566160068,
    1.57463669703976, 1.72177835457086, 1.79766772256782, 1.59000907896393,
    1.882096784056, 2.39819116843323, 2.36305207251141, 2.59145152780345,
    0.021066029607192, 0.021066029607192, -0.0297691100912388,
    0.0637691652631361, 0.411898608891399, 0.193912351274889, 0.387742231020071,
    0.387009323652292, 0.379445283918169, 0.374565459066037, 0.387009323652292,
    0.369684078147308, 0.362039122950075, 0.35821586731816, 0.341621973114356,
    0.347885141172337, 0.354555998679061, 0.00154517413206838,
    0.0194399400120143, 0.0194399400120143, 0.0194399400120143,
    0.0194399400120143, 0.019359024548905, 0.0196842424679406,
    -0.0601902120871469, 0.411329088516438, 0.61069078501842, 0.556763740989254,
    0.541878407910336, 0.601906789071264, 0.549524919174166, 0.579376500795113,
    0.622810987752237, 0.61630351723833, 0.617931162900106, 0.618254824752543,
    0.535370937396428, 0.480060550161413, 0.376761069036177, 0.635254852338492,
    0.651767831079186, 0.644854227183134, 0.6472941396092, 0.65014018541741,
    0.65022110088052, 0.647619357528235, 0.639811015338186, 0.637207715919303,
    0.636232062162197, 0.629480289192363, 0.621427644546387, 0.621508560009497,
    0.612235959150488, 0.513004036113094, 0.502348092048236, 0.539194193028344,
    0.534558670632139, 0.512107741752499, 0.524391331479133, 0.515280561546343,
    0.512026826289389, 0.522113249979286, 0.522438467898321, 0.492099837432121,
    0.458750218098294, 0.472090377045145, 0.470218428927442, 0.47379893817003,
    0.466315813899016, 0.463549127487316, 0.478759678485271, 0.467291467656123,
    0.465827208987164, 0.464524781244423, 0.460295392230363, 0.421253681280116,
    0.506983614444441, 0.613944520275374, 0.606948444849614, 0.617198255532327,
    0.677388467619474, -0.00585703467583622, -0.198547873671094,
    -0.199767829884127, -0.0155357689169896, -0.00431186054376785,
    -0.00870308048404676, -0.0350566243921133, -0.0974440025159979,
    -2.12050019150823, -1.90552181062615, -2.11635171795728, -1.75203606957395,
    -1.94480938009892, -1.64540193772865, -1.59733037230873, -1.71397001238081,
    -1.63710499062675, -1.47629795622976, -1.65296597746278, -1.69762042063294,
    -1.69932898175782, -1.59391325005896, -1.628319438613, -1.56999961857696,
    -1.59879463097769, -1.69005638089881, -1.6558944948007, -1.59367050366963,
    -1.55169871931486, -1.68175943379691, -1.60367445582982, -1.65931006098387,
    -1.67809956515781, -1.69981758666967, -2.1322127047933, -1.67224253048198,
    -2.22176744965614, -1.75472028445595, -1.98458399841695, -2.14880659899711,
    -2.37159332206924, -1.77058127129197, -2.123917313758, -1.7312952578858,
    -1.65613724119002, -1.68956777598696, -2.18443274376418, -1.73324656540001,
    -1.79034642922302, -1.75398893315477, -1.59171764008882, -1.67199822802605,
    -1.80425610854483, -1.87282418319699, -1.94407647273114, -2.38086592292825,
    -2.13026139727909, -2.03704678377715, -1.71323710501303, -1.90381480556786,
    -2.0411952573281, -1.83914934594413, -2.15002655521014, -2.15856780476797,
    -1.90357050311193, -1.84891055171499, -2.02265005561009, -1.7791225208498,
    -2.15075790651132, -1.74251916625902, -2.18467704622011, -1.91357445527212,
    -1.85671889390504, -1.5336421225087, -1.78937077546591, -1.86696714852115,
    -1.84110220952494, -1.7395906489211, -1.73812794631874, -2.20615076527604,
    -2.27349888371517, -1.94871355119394, -1.74276346871495, -1.84281077064983,
    -1.77814531102609, -1.94529642894417, -2.07316153345608, -1.65443023613174,
    -1.79059073167895, -2.05632489286295, -1.75496458691187, -2.0917067351741,
    -1.87868121787283, -1.89576216092188, -1.89478650716478, -1.99483225303306,
    -2.3308430499941, -2.15368642384924, -2.18052857266916, -1.69859607439004,
    -1.89600646337781, -2.13026139727909, -2.70077143059734, -2.19809812063007,
    -1.95432628341385, -2.27032606392132, -2.11366750307529, -3.53384280530574,
    -4.0575012294171, -5.02722170458457, -6.60332266029328, -5.80002039567612,
    -4.71805305644344, -5.83247527671595, -6.29317680232845, -5.98888953510606,
    -7.20531036482766, -6.49082993770555, -4.80272642038755, -4.97646592428265,
    -4.21586835138802, -4.7075604993714, -5.53819196165374, -6.20947909214145,
    -4.45500311014143, -3.18441093393474, -5.48768048380775, -5.77610676419412,
    -6.88369629603368, -5.56234833952507, -4.91497328445276, -6.80487996676541,
    -5.73194092593581, -7.55669354371287, -4.19488168117734, -5.61822824713505,
    -5.72437533013509, -5.71363847060712, -6.12700289423408, -5.79953334683086,
    -5.34249013379086, -5.55942137825375, -6.53206881469254, -3.23663097290562,
    -7.21799853186984, -6.30659787673841, -5.32418923452877, -4.35788588161107,
    -6.01890294765322, -5.04430264763362, -5.96375594741102, -5.66922832989289,
    -6.32343451733155, -3.84203424362316, -6.73728754587035, -4.55285324603957,
    -5.95399474164016, -6.49668697238139, -5.23780730945977, -7.02132260631645,
    -7.02351821628659, -5.67728097453887, -3.67512742816101, -6.35296088103346,
    1.07618499575356, 1.07456279632487, 1.06960516814283, 20.140600233644,
    28.7040322988636, 25.0940566011383, 22.5749971816522, 22.2970175543732,
    20.0552087449648, 19.9841525198221, 20.171718453476, 19.0907680037652,
    19.4160419411983, 19.057055820914, 18.7462571930112, 18.4892689041427,
    20.574260433912, 9.03006376063802, 1.0072559136506, 1.00579476711484,
    1.00697582166292, 1.00331284089062, 0.99814825585101, 1.01518796313521,
    1.72332664083613, 32.2231539361165, 32.9687486929041, 30.3733011860763,
    28.9705771707453, 30.4568036099052, 31.4257841754052, 29.9420093191708,
    30.9510142516779, 30.5332951756757, 30.4779489989089, 30.4824647041771,
    29.8279317427528, 39.290634145937, 34.8788577878619, 32.382614972906,
    32.3103076702183, 31.8550873891981, 31.9036934634945, 31.979298071304,
    32.0088983481693, 32.0073609543703, 31.5933460949052, 31.4211556553087,
    31.3600668147611, 31.1641167942758, 31.2887709553622, 30.9284489519033,
    30.3535391402784, 29.8571328885356, 26.4925994655198, 26.6007281993327,
    26.5787363100995, 26.60935658862, 27.7533814198577, 28.4241457119684,
    28.8429991605207, 27.4070811325157, 26.944678048078, 25.8719000601414,
    24.282620776412, 23.6271129392639, 23.503516125431, 23.4036781144536,
    23.1475378818812, 23.5620499046243, 23.735374716771, 23.2907745902839,
    23.1774275871338, 23.0942604956568, 23.2069088249043, 23.4951826107641,
    29.3715825354156, 30.499951002575, 30.471358278832, 33.521357736076,
    35.2905144283972, 24.8375164594501, 0.459848801116663, 1.47357639574941,
    1.74817157817718, 2.45634915754305, 2.18438761783283, 2.39452740962764,
    2.34788431334491, -48.3516823046795, -47.5907213902342, -47.9870118731118,
    -47.0776683380233, -48.0672862367827, -46.4669044178747, -47.2971725386064,
    -47.4750612940857, -47.0587683531209, -46.9435820792516, -47.2498338805542,
    -47.0835347090987, -47.5923015758648, -47.6747871101717, -46.7848547278631,
    -47.5147083149427, -46.3070745933015, -47.0460794080454, -47.6613660357618,
    -45.2491453708786, -45.9959507474796, -47.7487142781883, -47.8153364915593,
    -47.7639730672508, -47.4266256090818, -48.3624176081409, -47.1350233967684,
    -47.7625064744819, -47.4408955178211, -47.6296362817567, -46.6419899193773,
    -47.8077732298585, -48.2862823816214, -47.1646672434985, -47.1452693172846,
    -47.788503679139, -47.6533095009493, -47.4260187431085, -47.3410901842423,
    -47.2350636964036, -47.7066296790051, -47.6309807232976, -47.7278590956051,
    -47.3517087827088, -47.5531439379529, -47.8725538384105, -47.6712439465275,
    -47.2278645542868, -48.2404110943714, -47.3487755971711, -47.7399357284742,
    -47.5133630953685, -48.9065818778841, -47.1708891757916, -47.2517867441351,
    -47.1424606170748, -47.4289355899469, -47.2462860487102, -46.7233659781997,
    -46.7841249326285, -48.24249077738, -46.717628760652, -46.5870794412581,
    -48.0736264301372, -48.3810872951868, -46.98372704142, -47.5555830723457,
    -47.9258802407327, -46.972755993869, -47.3208364213994, -48.0122559415353,
    -46.9366404661567, -48.073503500876, -47.4685849449038, -47.0919538074286,
    -47.4042906071637, -47.7346910060048, -47.068774639381, -46.4863070122884,
    -47.444432457199, -46.699939395563, -47.9677220935265, -46.7487469804839,
    -48.6872076088619, -46.286210852352, -48.2602976254971, -46.5398652685338,
    -48.2799460784333, -47.1584507574385, -47.8153364915593, -46.6118504654357,
    -47.6825923402286, -46.4844786340355, -47.501168201438, -46.7971803313879,
    -48.2481003974667, -46.8530617950645, -46.8589141615405, -48.262984174479,
    -44.0388196540224, -53.7366382739701, -47.803256746557, -48.9733161280502,
    -48.2598066864854, -45.3763981631841, -49.2745519486734, -53.1109812404097,
    -50.5672324902424, -48.2013709394878, -49.5177169199163, -51.3406372692739,
    -48.5567897770956, -47.3916039961882, -52.348671994257, -47.8604818737412,
    -48.747241436256, -41.9440630283157, -43.5102860732587, -53.4773680134353,
    -45.4126311739245, -51.7989012165603, -44.2978502803011, -44.8993454897572,
    -50.2104847717596, -45.6019889182663, -51.7888964863668, -44.0498015940396,
    -47.4527325164341, -42.3601085546915, -49.4021696385962, -49.3661809303117,
    -48.0045759748396, -48.7895757841281, -47.8358384470246, -48.832887341824,
    -41.5532709067631, -46.6473529029082, -42.849850173194, -42.1198806591528,
    -45.7665834187635, -48.7892108865108, -45.4793724263903, -46.7112854551642,
    -49.5079596043119, -49.8553187949741, -41.8591438058491, -49.3843612344125,
    -42.7527274984305, -45.964599117658, -48.6158323900665, -51.22449479248,
    -47.8437658283094, -52.970546229917, -49.7462439726694, -44.9276511192127,
    -52.0071737284372, 0.0, 0.0, 0.0, 12.0, 13.0, 15.0, 13.0, 15.0, 17.0, 16.0,
    13.0, 17.0, 16.0, 17.0, 15.0, 18.0, 15.0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 12.0, 13.0, 12.0, 15.0, 11.0, 16.0, 17.0, 16.0, 15.0, 16.0, 17.0, 17.0,
    12.0, 17.0, 15.0, 18.0, 17.0, 19.0, 21.0, 20.0, 19.0, 19.0, 19.0, 17.0, 15.0,
    16.0, 15.0, 17.0, 15.0, 15.0, 20.0, 19.0, 18.0, 17.0, 15.0, 13.0, 18.0, 14.0,
    15.0, 14.0, 17.0, 17.0, 16.0, 17.0, 13.0, 19.0, 17.0, 20.0, 20.0, 17.0, 15.0,
    15.0, 18.0, 14.0, 12.0, 13.0, 15.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.0, 7.0, 8.0, 8.0, 8.0, 8.0, 8.0,
    7.0, 7.0, 8.0, 8.0, 7.0, 8.0, 8.0, 7.0, 8.0, 8.0, 8.0, 7.0, 7.0, 8.0, 8.0,
    8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0,
    8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0, 8.0,
    8.0, 8.0, 8.0, 8.0, 10.0, 0.0, 0.0, 0.0, 46.0, 42.0, 43.0, 43.0, 43.0, 43.0,
    44.0, 45.0, 44.0, 44.0, 47.0, 44.0, 47.0, 46.0, 20.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 43.0, 43.0, 41.0, 45.0, 41.0, 43.0, 43.0, 47.0, 43.0, 43.0,
    47.0, 44.0, 42.0, 46.0, 46.0, 46.0, 44.0, 46.0, 47.0, 44.0, 46.0, 46.0, 44.0,
    44.0, 45.0, 47.0, 45.0, 44.0, 46.0, 45.0, 47.0, 45.0, 42.0, 43.0, 44.0, 44.0,
    46.0, 43.0, 44.0, 44.0, 46.0, 44.0, 40.0, 45.0, 44.0, 44.0, 45.0, 44.0, 45.0,
    41.0, 45.0, 45.0, 45.0, 47.0, 39.0, 45.0, 46.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 39.0, 40.0, 44.0, 43.0,
    43.0, 44.0, 44.0, 39.0, 39.0, 46.0, 44.0, 42.0, 47.0, 44.0, 38.0, 45.0, 45.0,
    45.0, 39.0, 39.0, 45.0, 44.0, 44.0, 43.0, 44.0, 44.0, 43.0, 42.0, 43.0, 44.0,
    43.0, 43.0, 45.0, 43.0, 42.0, 42.0, 42.0, 42.0, 42.0, 42.0, 42.0, 44.0, 42.0,
    44.0, 42.0, 42.0, 42.0, 42.0, 42.0, 41.0, 41.0, 41.0, 42.0, 41.0, 43.0, 38.0,
    42.0, 0.0, 0.0, 0.0, 3.83333333333333, 3.23076923076923, 2.86666666666667,
    3.30769230769231, 2.86666666666667, 2.52941176470588, 2.75, 3.46153846153846,
    2.58823529411765, 2.75, 2.76470588235294, 2.93333333333333, 2.61111111111111,
    3.06666666666667, 3.33333333333333, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    3.58333333333333, 3.30769230769231, 3.41666666666667, 3.0, 3.72727272727273,
    2.6875, 2.52941176470588, 2.9375, 2.86666666666667, 2.6875, 2.76470588235294,
    2.58823529411765, 3.5, 2.70588235294118, 3.06666666666667, 2.55555555555556,
    2.58823529411765, 2.42105263157895, 2.23809523809524, 2.2, 2.42105263157895,
    2.42105263157895, 2.31578947368421, 2.58823529411765, 3.0, 2.9375, 3.0,
    2.58823529411765, 3.06666666666667, 3.0, 2.35, 2.36842105263158,
    2.33333333333333, 2.52941176470588, 2.93333333333333, 3.38461538461538,
    2.55555555555556, 3.07142857142857, 2.93333333333333, 3.14285714285714,
    2.70588235294118, 2.58823529411765, 2.5, 2.64705882352941, 3.38461538461538,
    2.31578947368421, 2.64705882352941, 2.2, 2.25, 2.41176470588235, 3.0, 3.0,
    2.5, 3.35714285714286, 3.25, 3.46153846153846, 3.06666666666667, 2.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    5.57142857142857, 5.71428571428571, 5.5, 5.375, 5.375, 5.5, 5.5,
    5.57142857142857, 5.57142857142857, 5.75, 5.5, 6.0, 5.875, 5.5,
    5.42857142857143, 5.625, 5.625, 5.625, 5.57142857142857, 5.57142857142857,
    5.625, 5.5, 5.5, 5.375, 5.5, 5.5, 5.375, 5.25, 5.375, 5.5, 5.375, 5.375,
    5.625, 5.375, 5.25, 5.25, 5.25, 5.25, 5.25, 5.25, 5.25, 5.5, 5.25, 5.5, 5.25,
    5.25, 5.25, 5.25, 5.25, 5.125, 5.125, 5.125, 5.25, 5.125, 5.375, 4.75, 4.2,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 2.0, 4.0, 0.0, 5.0, 1.0,
    2.0, 2.0, 2.0, 2.0, 1.0, 0.0, 1.0, 1.0, 3.0, 1.0, 2.0, 1.0, 0.0, 1.0, 2.0,
    1.0, 1.0, 3.0, 4.0, 2.0, 4.0, 2.0, 1.0, 4.0, 4.0, 1.0, 4.0, 2.0, 3.0, 6.0,
    5.0, 3.0, 6.0, 6.0, 1.0, 4.0, 4.0, 4.0, 5.0, 4.0, 3.0, 4.0, 2.0, 2.0, 4.0,
    4.0, 5.0, 5.0, 3.0, 5.0, 0.0, 0.0, 0.0, 31.0, 50.0, 47.0, 46.0, 48.0, 20.0,
    32.0, 37.0, 10.0, 17.0, 16.0, 5.0, 2.0, 28.0, 9.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 49.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 44.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.0, 9.0, 9.0, 7.0, 12.0, 5.0,
    11.0, 8.0, 8.0, 9.0, 7.0, 4.0, 10.0, 9.0, 8.0, 8.0, 10.0, 8.0, 8.0, 11.0,
    9.0, 8.0, 11.0, 10.0, 8.0, 7.0, 9.0, 8.0, 8.0, 10.0, 6.0, 7.0, 10.0, 7.0,
    9.0, 6.0, 5.0, 5.0, 8.0, 4.0, 5.0, 12.0, 7.0, 5.0, 6.0, 5.0, 6.0, 8.0, 8.0,
    7.0, 11.0, 7.0, 6.0, 5.0, 5.0, 5.0, 6.0, 50.0, 50.0, 50.0, 19.0, 0.0, 3.0,
    4.0, 2.0, 30.0, 18.0, 13.0, 40.0, 33.0, 34.0, 45.0, 48.0, 22.0, 41.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 49.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 6.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 39.0, 41.0, 39.0, 39.0, 38.0, 40.0, 38.0,
    40.0, 40.0, 39.0, 41.0, 45.0, 40.0, 40.0, 41.0, 39.0, 39.0, 40.0, 41.0, 39.0,
    40.0, 40.0, 38.0, 39.0, 39.0, 39.0, 39.0, 38.0, 40.0, 39.0, 40.0, 39.0, 39.0,
    39.0, 39.0, 41.0, 39.0, 40.0, 39.0, 40.0, 39.0, 37.0, 39.0, 41.0, 40.0, 40.0,
    40.0, 39.0, 38.0, 41.0, 37.0, 39.0, 40.0, 40.0, 40.0, 42.0, 39.0,
    0.000813822830887934, 1.55606659835161E-6, 0.170404853175033,
    0.0579930460504091, 0.0258633829296159, 0.144295611722893,
    0.0528704748089499, 0.0773520705993142, 0.0084603340947188,
    0.0156182404466972, 0.00658838597701661, 0.0055318167568007,
    0.0327785428922662, 0.00821603163879259, 0.0176504634240198,
    0.127376499600054, 0.0310715378339792, 0.313236206228967, 0.0,
    0.000731351301180405, 0.000325217919035503, 0.000245858522524567,
    0.00105656922021591, 0.000569520374961718, 0.0566937304408655,
    0.0796317081657595, 0.0573441662789365, 0.0505130339125918,
    0.0230982525845146, 0.000894738293997333, 0.034161886098116,
    0.173902890887913, 0.00227652543324874, 0.00114059681652179,
    0.00715790635197833, 0.00129931560954399, 0.0138272077921043,
    0.0463630042950426, 0.176096944791454, 0.0784101958861284,
    0.00569364768301917, 0.00203377904392088, 0.0030094328010275,
    0.00121995621303284, 0.000977209823704861, 0.00268421488199189,
    0.00105812528681415, 0.00121995621303284, 0.00544934522709306,
    8.24715297075285E-5, 0.00170856112488527, 0.0029285173379181,
    0.00488138091872958, 0.0288743717972415, 0.0109002465207845,
    0.00105656922021591, 0.00626472412457935, 0.0259474105259215,
    0.0579136866538981, 0.0963874332957819, 0.0466882222140782,
    0.0115506823588556, 0.0393684849358811, 0.0511619136840645,
    0.0117934287481836, 0.010411641608932, 0.0107353034613693,
    0.00300943280102739, 0.000406133382144902, 0.0145585590932847,
    0.000242746389327975, 0.00740065274130641, 0.0022780814998471,
    0.0335939217897524, 0.0239945469451102, 0.0280621050329518,
    0.225063248505379, 0.00178947658799467, 0.0259474105259216,
    0.166500682080008, 0.0361956651420367, 0.266709814941107, 0.0989082611849569,
    0.121765323446742, 0.00748312427101383, 0.0122011181969267,
    0.00227808149984698, 0.00276668641169942, 0.258738085758241,
    0.172765406204588, 0.0422145307440918, 0.0668610695938711, 0.343330534239241,
    0.0876018812820276, 0.289403490210075, 0.0788178853348715, 0.22010095212354,
    0.376028161668398, 0.0412388769869849, 0.0571014198896084,
    0.0456316529938622, 0.164711205492014, 0.488033835410877, 0.296965973877601,
    0.0700323333211173, 0.134209188032997, 0.0714965919900761,
    0.0790606317241995, 0.0836961541204047, 0.0822334515180438,
    0.025864938996214, 0.16495550794794, 0.0490487752436326, 0.206438687390851,
    0.256950165236845, 0.0300134125471651, 0.132500626908111, 0.125179333563316,
    0.0468516092068951, 0.0453857944713377, 0.0714981480566743,
    0.346503354033086, 0.189600490731121, 0.185697875702695, 0.16007412702921,
    0.166173908094375, 0.285743621570977, 0.0663724646820187, 0.145434652472816,
    0.202777262685155, 0.123715074894357, 0.368462565867677, 0.192774866591564,
    0.560994686069913, 0.206927292302704, 0.0463630042950425, 0.0858933201571422,
    0.00292696127131986, 0.505357524849263, 0.119323854954078, 0.061004034918035,
    0.136162051613808, 0.0568571174336823, 0.0646654596237317, 0.127133753210726,
    0.287937675474518, 0.0927260085900852, 0.185453573246769, 0.294039012606281,
    0.359435823531193, 0.211561258632311, 0.357727262406307, 0.399454744305145,
    0.121030860012365, 0.232059323931139, 0.175935113865236, 0.0175695479609105,
    0.394087870607759, 0.0329419298850832, 0.427759595727424, 0.26792977115414,
    0.0549042538528707, 0.061248337373961, 0.069788030865191, 0.254752999200107,
    0.561480178848569, 0.196921784075916, 0.264758507426894, 0.501940402599493,
    0.15299869220694, 0.120543811167111, 0.284279362902018, 0.121278274601488,
    0.0597840787050021, 0.053684297639838, 0.428492503095203, 0.254752999200107,
    0.160075683095809, 0.0751564606291748, 0.328200898704397, 0.3208811614262,
    0.346259051577159, 0.110538302940323, 0.288183533997043, 0.197410388987769,
    0.015859430769427, 0.208877043750319, 0.618581598738177, 0.548790455739789,
    1.48386355203049, 0.455090349459196, 0.215465429727336, 1.19324010560738,
    1.0582995662732, 0.06759397696165, 0.9709404313805, 0.763039041387288,
    0.162026990610022, 1.11320070799287, 1.17859751891779, 1.17054487427181,
    0.31868555145606, 0.32722524494729, 1.07049912840353, 0.869430426843259,
    0.788172629082323, 0.801349401036355, 0.415071428685244, 1.14492268166492,
    0.978504471114624, 0.136162051613808, 0.160318429485137, 0.45557739830445,
    0.689832332205728, 1.1685935667576, 0.334056377313634, 0.269638332279026,
    1.29304154901973, 1.21617652726567, 1.41260970642974, 0.428736805551129,
    0.65225487992444, 1.57512530195161, 0.83282707225248, 1.0968526723116,
    0.955323747000402, 0.453626090790237, 0.292819056393248, 0.384080806314374,
    1.40065289068874, 0.00390417109502494, 1.23008465052089, 0.714967475967361,
    0.882362896341367, 0.861131923674759, 0.430201064220089, 1.19299424708485,
    1.10783383429549, 0.427272546882171, 1.18738307093154, 1.72861104300381,
    0.25597295541314, 1.11051804917748, 1.16615365433153, 0.138601964039874,
    0.000813822830887934, 3.11213319670323E-7, -0.170404853175033,
    0.0193310153501364, 0.00287370921440176, -0.144295611722893,
    -0.0528704748089499, 0.0773520705993142, -0.0042301670473594,
    0.000600701555642202, -0.00658838597701661, 0.000691477094600088,
    0.0327785428922662, 0.00821603163879259, -0.0176504634240198,
    -0.00310674389268423, 0.0310715378339792, 0.313236206228967, 0.0,
    0.000146270260236081, 0.000162608959517752, -1.75613230374691E-5,
    -0.00105656922021591, 0.000569520374961718, 0.00515397549462413,
    -0.0113759583093942, 0.00819202375413378, 0.0101026067825184,
    -0.000745104922081116, -0.000127819756285333, -0.034161886098116,
    0.00445904848430546, 0.000455305086649749, 0.00114059681652179,
    0.00715790635197833, 0.000259863121908799, -0.0138272077921043,
    0.00331164316393161, -0.0092682602521818, -0.00326709149525535,
    -0.000813378240431309, 0.00203377904392088, -0.0030094328010275,
    0.00121995621303284, 7.51699864388355E-5, 0.00268421488199189,
    -3.20644026307318E-5, 0.00121995621303284, 0.00544934522709306,
    8.24715297075285E-5, 0.00170856112488527, -0.000266228848901646,
    -0.000406781743227465, -0.0288743717972415, 0.00545012326039224,
    0.000150938460030844, -0.00626472412457935, -0.00288304561399128,
    0.0144784216634745, 0.0963874332957819, 0.0466882222140782,
    -0.0115506823588556, 0.0393684849358811, -0.0511619136840645,
    -0.0117934287481836, -0.010411641608932, -0.0107353034613693,
    0.00300943280102739, -5.80190545921289E-5, -0.0145585590932847,
    -0.000242746389327975, 0.00740065274130641, -0.0022780814998471,
    0.0335939217897524, -0.000749829592034693, -0.000876940782279744,
    -0.225063248505379, -0.000127819756285333, 0.00518948210518433,
    0.00370001515733352, 0.00603261085700611, -0.266709814941107,
    -0.00494541305924785, 0.0608826617233709, -0.00748312427101383,
    0.0122011181969267, 0.00113904074992349, -0.000172917900731213,
    0.258738085758241, -0.010162670953211, 0.00211072653720459,
    0.00180705593496949, 0.11444351141308, 0.00461062533063303,
    -0.0180877181381297, 0.0197044713337179, -0.00489113226941201,
    -0.0235017601042749, 0.00294563407049892, 0.00439241691458526,
    -0.00285197831211639, 0.0149737459538194, -0.488033835410877,
    0.00802610740209731, 0.0700323333211173, 0.134209188032997,
    0.003249845090458, 0.0790606317241995, 0.00288607428001395,
    0.0822334515180438, -0.000957960703563482, 0.0126888852267646,
    0.00222948978380148, 0.206438687390851, -0.00734143329248128,
    0.00230872404208962, 0.0101923559160085, -0.0113799394148469,
    -0.00141974573354228, -0.00216122830815894, -0.00649983164151585,
    0.346503354033086, -0.0189600490731121, -0.0309496459504491,
    -0.16007412702921, -0.00664695632377499, 0.285743621570977,
    -0.0110620774470031, 0.00338220122029806, 0.00844905261188144,
    -0.123715074894357, 0.0153526069111532, 0.0642582888638546,
    0.0280497343034956, 0.00985368058584304, -0.0463630042950425,
    0.0026841662549107, -6.36295928547796E-5, 0.505357524849263,
    -0.119323854954078, 0.00508366957650291, 0.045387350537936,
    0.0142142793584206, 0.00170172262167715, -0.00385253797608259,
    0.287937675474518, -0.0154543347650142, 0.0618178577489229,
    -0.014701950630314, 0.359435823531193, 0.211561258632311, 0.357727262406307,
    0.399454744305145, 0.00605154300061825, -0.232059323931139,
    -0.00651611532834206, 0.0175695479609105, -0.0145958470595466,
    -0.00183010721583795, -0.0237644219848569, -0.00992332485756076,
    0.00499129580480643, 0.00157047018907592, 0.00303426221153005,
    0.0159220624500067, -0.561480178848569, -0.0123076115047448,
    -0.0126075479727092, 0.0125485100649873, 0.00402628137386684,
    0.120543811167111, 0.284279362902018, -0.00466454902313414,
    -0.0199280262350007, 0.00185118267723579, 0.0142830834365068,
    -0.0149854705411828, 0.160075683095809, 0.0044209682723044,
    -0.0193059352179057, -0.0133700483927583, 0.346259051577159,
    -0.00650225311413667, 0.144091766998521, -0.00564029682822197,
    -0.000932907692319234, -0.0261096304687898, -0.0134474260595256,
    -0.0274395227869895, 1.48386355203049, 0.455090349459196,
    0.00567014288756146, 1.19324010560738, 1.0582995662732, 0.0168984942404125,
    0.9709404313805, 0.0305215616554915, 0.0135022492175018, 1.11320070799287,
    1.17859751891779, 1.17054487427181, 0.31868555145606, 0.0181791802748495,
    1.07049912840353, 0.869430426843259, 0.788172629082323, -0.200337350259089,
    0.0122079831966248, 1.14492268166492, 0.978504471114624, -0.0151291168459786,
    -0.160318429485137, 0.09111547966089, 0.0176880085180956, -1.1685935667576,
    0.00927934381426762, -0.0112349305116261, 1.29304154901973,
    -1.21617652726567, 1.41260970642974, -0.428736805551129, -0.05017345230188,
    1.57512530195161, 0.83282707225248, -0.0997138793010541, 0.955323747000402,
    -0.0151208696930079, -0.00791402855116887, -0.0103805623328209,
    -1.40065289068874, 0.000108449197084026, 1.23008465052089,
    -0.0549974981513354, 0.882362896341367, 0.861131923674759,
    -0.0165461947776957, -1.19299424708485, 0.0307731620637636,
    0.427272546882171, 1.18738307093154, 1.72861104300381, -0.0196902273394723,
    1.11051804917748, 1.16615365433153, -0.0126001785490794, 49.6329708826827,
    49.6389088328216, 49.8800057915554, 43.6893115353881, 39.9309765580185,
    42.0766928087515, 43.9368226006399, 44.298703004797, 45.3398344882916,
    45.2909491000408, 45.061983236507, 45.5757902029848, 45.5023531959467,
    45.7769460442745, 45.802244575029, 45.873985469475, 44.975031791063,
    47.6817333914663, 49.6293141261768, 49.6349237462635, 49.6302182008703,
    49.6300517017443, 49.6505435427767, 49.6184185478558, 50.1111641529261,
    36.5489124887067, 36.6904336336849, 38.0975971070207, 38.790200793838,
    38.7724631906844, 37.0789943557357, 38.9247289755236, 38.2102516605347,
    38.4826520110256, 38.4941264461212, 38.470209702506, 38.9498641192852,
    10.4820427300369, 22.1111275044456, 37.2740473038271, 37.2848526302854,
    37.6820386294902, 37.6439803526301, 37.5551491787354, 37.5257006183634,
    37.5328569686488, 37.8831073310505, 38.0239080171937, 38.0751477342077,
    38.2444991302957, 38.212442602305, 38.5116913258823, 38.9849129633444,
    39.3253289848451, 41.5978428784036, 41.5086211317927, 41.5349731196342,
    41.5109754605559, 40.5119729241424, 39.9335004980409, 39.5774709042933,
    40.7599725943061, 41.1292520951379, 42.0401610332242, 43.0783655554476,
    43.4722931511957, 43.5093897788981, 43.5503065500993, 43.6716610719641,
    43.4436615257878, 43.3533054066269, 43.6013004085907, 43.6394442691138,
    43.6359540117339, 43.5472099775687, 43.5930050175554, 39.6601929607232,
    38.9363636854788, 38.9720691896424, 36.2564933414292, 34.4709520458616,
    39.9614739072778, 27.3809773372967, 32.8380717760491, 33.317561249769,
    35.8896630894231, 39.3181679663599, 39.9526930234638, 38.6939658550688,
    1.57145765297952, 3.01310667426592, 2.47920155160504, 2.54142554273611,
    3.50944056891189, 2.83083214507938, 2.84181486312988, 1.71128579749882,
    4.01382244000405, 2.45382210538749, 3.73393274097571, 2.84449596587867,
    3.40109943200831, 2.17027720591984, 3.74271829298947, 4.4620801010303,
    2.9247765538159, 3.13585231966956, 3.23346437737817, 3.75980079210512,
    3.80299408873952, 2.38355169387684, 2.66855929594684, 4.1043512825574,
    2.67710210157126, 2.69832996210467, 3.97746338786919, 3.05801164415839,
    3.46259051577159, 1.8315884183432, 2.92795404180954, 3.33204275244429,
    2.36085490647467, 3.77810324743381, 3.13755465652805, 3.86741213377413,
    2.46310404264608, 1.93676140365274, 3.70416674301767, 3.29788086634618,
    3.64632930362709, 3.32203568815091, 3.70758386526744, 3.98649324233887,
    3.65218789436952, 3.39401777291964, 3.21272512175661, 4.27125965408615,
    3.29032305087845, 2.9247781098825, 2.78276057965946, 3.28885412400969,
    3.84081117527693, 2.01679146486765, 3.92671071970047, 2.27446986927248,
    3.25444793545565, 3.9711208604147, 3.32960128395163, 3.29641349554402,
    2.35817535979247, 1.56195008606418, 3.76102852865115, 4.26832646854843,
    2.45480087127779, 2.50677660779274, 4.69633036673178, 2.08023541221175,
    3.72173784704518, 3.06777440599585, 3.03751357885956, 3.65341251878235,
    3.30861461374095, 2.82253831011068, 2.67343134046598, 3.95257721476328,
    3.88327312061015, 3.15415633106485, 3.38596512827367, 3.74760434210799,
    3.36937434620306, 2.91257543561898, 2.68832445387789, 3.32349527862007,
    5.01233637150565, 2.90257615165858, 3.88961876019784, 4.0179662453552,
    4.44207375277652, 3.82007503178857, 3.57361897969037, 3.25102614500609,
    4.16267577079324, 3.04239495977829, 3.93671622792726, 4.55284857783978,
    3.99185389176987, 3.41159198908035, 3.50602500272872, 1.77692846694626,
    8.39122606495704, 5.83808800893586, 7.08696060756409, 7.28974565058224,
    3.01018282512779, 16.2622108280525, 4.66168609998821, 8.77872554958779,
    7.56377987100133, 7.85512377625944, 4.02138647973817, 7.02766668983754,
    6.56915921812857, 9.14158316356855, 7.60549334830079, 8.95051530203527,
    3.35083070055164, 2.29692235421872, 10.498554152711, 9.79505021937302,
    3.1021774824167, 9.09301676897037, 9.64913474232135, 9.9631863274478,
    5.72511290570267, 7.23362299658293, 2.96235089396399, 5.552591801954,
    5.56479136408433, 9.62717241835355, 5.61407977358409, 9.37534482435818,
    10.0664033370433, 6.77364504193861, 10.0951952373109, -3.58606128821003,
    1.75716330901521, 8.55861681713125, 6.49595250894701, 2.39404269488228,
    11.8870059252028, 5.78074851085672, 15.0396917733323, 10.8923899405299,
    -0.336491621539905, -1.20690392640667, 18.0406035536677, -1.23618131945285,
    12.872344640817, 7.14821361313784, 8.86583726989136, 8.36878447247699,
    17.6291888814884, 20.6408468577513, 3.13390412428854, 24.5914649846674,
    0.000457586710872458, 0.000303135455371289, 0.0488618929163801,
    0.106563358873831, 0.09539061778589, 0.0758998769453926, 0.0367837609197251,
    0.0269824434364968, 0.0207068423307799, 0.0154282031596201,
    0.0158153056500143, 0.00869195280310732, 0.0338735220937851,
    0.0170885792806723, 0.0316075762429317, 0.0272728372914003,
    0.0414353995264101, 0.0907028604241813, 0.000309609277930563,
    0.000380195687508055, 0.000334606167738764, 0.000295933857113724,
    0.00253860175991257, 0.000506207293896572, 0.035849211285308,
    0.0792567856606008, 0.0279074107487349, 0.0547612949675654,
    0.0276560005756694, 0.016707439261804, 0.0875915180858257, 0.082984148937124,
    0.0113832573938204, 0.0068662333370892, 0.00602100471424646,
    0.0042023952666758, 0.0307999417667166, 0.5098853791597, 0.50118752211284,
    0.0300298347980431, 0.00982486058311123, 0.00750765011616282,
    0.00520476026238853, 0.00469764169439378, 0.0056004069603083,
    0.00760595605133476, 0.00475316135061478, 0.00588506376484696,
    0.0071530134000876, 0.00583519453141601, 0.0163014756955883,
    0.0067476667631236, 0.00590187192853539, 0.03405837920322,
    0.0188702304015852, 0.00343322515120996, 0.00669087482866091,
    0.00781967886317348, 0.0431202195082375, 0.016484697605125,
    0.0517561618571037, 0.0146068722166235, 0.0197918605754512,
    0.0244363535411138, 0.011379349282046, 0.00638190154099958,
    0.00706584309737915, 0.00486251314001403, 0.00557030055039891,
    0.0185265399795573, 0.00488533344427038, 0.0111637743173804,
    0.00294754366641072, 0.0282505877222958, 0.0203697635996663,
    0.0280017096469253, 0.0561875016553133, 0.00981589902559519,
    0.00840727502857113, 0.0638869967637412, 0.0661619407098012,
    0.201487132756789, 0.154736151205412, 0.123758145109271, 0.0366794785368608,
    0.0241248147704377, 0.0260979057285628, 0.00854114055182117,
    0.0824291950978087, 0.248977205532711, 0.279122733783075, 0.261119981098315,
    0.243921387018905, 0.243423727477827, 0.242721789870969, 0.191423190153591,
    0.214062978563675, 0.21568038427691, 0.216036420887172, 0.212936442288901,
    0.20405900723529, 0.209440125709208, 0.22741133612545, 0.228967085809854,
    0.198278449807156, 0.195610384013653, 0.196807468118312, 0.193762813832384,
    0.15245634247781, 0.173937874442051, 0.212665543665674, 0.246100626541833,
    0.191872686187934, 0.256509309268983, 0.264016845195156, 0.252697734893395,
    0.227803094899961, 0.234503123427523, 0.241984089276819, 0.243459827956075,
    0.276858574167823, 0.263243022511512, 0.255839338910219, 0.249905024069041,
    0.222268843548647, 0.24573863834399, 0.239878748532881, 0.237277831817301,
    0.234073565806148, 0.218467177940863, 0.238858002481157, 0.253368614903099,
    0.20519759268978, 0.267162863887694, 0.260370085838426, 0.273186499043378,
    0.24173827514082, 0.235063188357818, 0.266616162453232, 0.245438019555577,
    0.258779953940133, 0.258869737383332, 0.243640485303786, 0.250023706649175,
    0.272239348694113, 0.275892840202083, 0.255530497826556, 0.23355016317897,
    0.256455269493767, 0.281144368499814, 0.279633047245727, 0.267993039423829,
    0.219832631649677, 0.246541008538637, 0.254121116086334, 0.255314186915077,
    0.266608828973695, 0.245525789451013, 0.268651004504778, 0.258841970608493,
    0.255578984280611, 0.257934940927949, 0.263777722711468, 0.228789853350624,
    0.262722813456036, 0.271799868657197, 0.266550656038193, 0.247138877835497,
    0.276040656473104, 0.266617254174192, 0.266342640007437, 0.265216798393025,
    0.264073750855391, 0.280505749024655, 0.264571842728576, 0.278227217766253,
    0.262743442324155, 0.248034504043917, 0.296883283359563, 0.275104315151692,
    0.261949059193746, 0.274270810852804, 0.285719194098502, 0.248127773925314,
    0.254989348663747, 0.270916427490574, 0.249464633203871, 0.316275850640018,
    0.603435815903032, 0.580033214269639, 0.692103741490943, 0.710573062938783,
    0.613439123909371, 0.719599569845222, 0.457036239278408, 0.81326529921808,
    0.758370047512972, 0.559135306265312, 0.747980279754979, 0.775806299038652,
    0.606427730537019, 0.650303313425062, 0.572870906389807, 0.79666513915646,
    0.672571665139563, 0.748206900368058, 0.727094673290473, 0.641552666097892,
    0.677470662492181, 0.81549596458439, 0.468975765023907, 0.585506066546311,
    0.628625167860333, 0.755487727718301, 0.628319034050446, 0.781137000279791,
    0.765508680370343, 0.732851017482296, 0.671212964298161, 0.828518971847252,
    0.570186812099424, 0.674840320504718, 0.765354916012668, 0.676322973783003,
    0.62249120940858, 0.840866822963961, 0.665849614270078, 0.547306563565059,
    0.695385147293785, 0.683726983566316, 0.554866197849651, 0.680725428583031,
    0.918285033276039, 0.956028946021618, 1.01737957081946, 0.751771359158652,
    0.849332986632414, 0.866265074864601, 0.651492623550378, 0.747511531181683,
    0.872355565772993, 0.651150018510031, 1.08274458055682, 1.26709877727257,
    1.46587945999619, 0.993389804193542, 0.993308888730433, 1.21056379504576,
    1.2108905690314, 1.01535212816133, 1.06065700716956, 0.991030807230586,
    0.964352045403484, 0.955242831537292, 0.946457279523538, 0.963782525028522,
    0.931816248900546, 0.999083451876562, 0.951175273449451, 1.00542753539765,
    0.964840650315336, 0.98297816258461, 1.32704160419163, 0.993226417200725,
    0.993797493642285, 0.993389804193542, 0.993227973267324, 1.00892557311053,
    0.993064586274507, 1.16607273886842, 1.02820368219633, 0.827052509106351,
    0.926284432143745, 0.875611123371533, 0.811516740189361, 1.02047625546939,
    0.973868948718419, 0.784187542524188, 0.786383152494328, 0.785082280818185,
    0.776702862186576, 0.890659843443268, 1.07236952045463, 1.1621685677734,
    0.819650300298446, 0.767430261327567, 0.774995857128289, 0.763852864218176,
    0.761086177806477, 0.761655698181439, 0.770683996584521, 0.767837950776311,
    0.773937731841475, 0.780200899899456, 0.77702963617221, 0.814932306372533,
    0.782559896862412, 0.792808151478527, 0.89757344733932, 0.878295338253525,
    0.836325109965359, 0.853325137551308, 0.862922956329352, 0.924658342548567,
    0.897329144883394, 0.891065976825413, 0.879270992010631, 0.893586804714588,
    0.892285933038446, 0.896026717140654, 0.893750191707405, 0.889439887230235,
    0.879922983915301, 0.883908070473435, 0.916280479983556, 0.880492504290262,
    0.901720364823673, 0.878703027702268, 0.965735388609334, 0.902697574647378,
    0.941577454671407, 0.994040240031613, 0.804114531381457, 0.810947219814399,
    0.911481570594534, 0.879922983915301, 1.1971427206358, 0.802894575168424,
    0.933117120576688, 0.771253516959483, 0.756287268417455, 0.833152290171515,
    0.826563904194498, 0.896190104133471, 0.40897164762008, 0.641762322852399,
    0.432396674190228, 0.433616630403261, 0.50413756863623, 0.424832634456105,
    0.439473665079097, 0.396526226967227, 0.52609833653742, 0.42727254688217,
    0.564165949797159, 0.349918920216258, 0.463142994105172, 0.688369629603368,
    0.571485687075356, 0.55147622668838, 0.379200981462243, 0.361875735957258,
    0.42092846336108, 0.439960713924351, 0.39555057321012, 0.501453353754239,
    0.39213345096035, 0.452406134577204, 0.456310305672229, 0.566361559767298,
    0.458994520554221, 0.46045877922318, 0.585395366397168, 0.444597792387155,
    0.493400709108263, 0.558308915121322, 0.497060577747361, 0.506577481062296,
    0.503161914879124, 0.468267121413229, 0.472659897420106, 0.439229362623171,
    0.516094384377231, 0.410922955134293, 0.481443893367262, 0.680561287413318,
    0.63468533196353, 0.350651827584037, 0.606867529386505, 0.61492017403248,
    0.642250927764252, 0.610284651636275, 0.521707116597141, 0.580758287934364,
    0.595399318557357, 0.642495230220178, 0.624437077347415, 0.47802832718409,
    0.477295419816311, 0.567337213524405, 0.494620665321296, 0.505357524849263,
    0.476808370971057, 0.532199673669182, 0.516825735678411, 0.446793402357294,
    0.541227972072265, 0.508773091032436, 0.437276499042359, 0.544887840711363,
    0.638102454213301, 0.489252235557312, 0.644690840190317, 0.444842094843081,
    0.663967393209515, 0.439229362623171, 0.549036314262314, 0.732781326384197,
    0.464118647862278, 0.534883888551174, 0.900419493147531, 0.562213086216347,
    0.508530344643108, 0.639811015338186, 0.569777125950471, 0.489740840469164,
    0.595643621013283, 0.607111831842431, 0.60857609051139, 0.549769221630093,
    0.613700217819447, 0.534395283639322, 0.605891875629398, 0.615653081400259,
    0.606867529386505, 0.49779348511514, 0.538788059646199, 0.595887923469209,
    0.60638048054125, 0.552209134056158, 0.558797520033175, 0.482419547124369,
    0.578074073052372, 1.36673530704655, 1.62124400379073, 2.27447453747228,
    1.63246791216395, 1.46141262315085, 1.45311567604894, 1.52314956543666,
    1.56707265730564, 1.79132208298013, 1.84671494174485, 1.82158135404982,
    1.38454760139679, 1.43749899166884, 1.28889307546878, 1.59196194254474,
    1.73617508273793, 1.38503620630864, 1.55560289040989, 1.48361769350796,
    1.44067181146269, 1.76155297288889, 1.36478244346574, 1.45043146116695,
    1.75276897694173, 1.92577557346905, 1.71031014374171, 1.63734929308268,
    1.54267042091178, 2.00093203409822, 1.41163405267263, 1.95481333225911,
    2.11635171795728, 1.86160027482377, 1.62368391621679, 1.80547606475787,
    1.63344512198765, 1.28498890437376, 1.3984572807186, 1.70811297770498,
    1.63808064438386, 2.03387552004991, 1.89185798982686, 1.34086881198373,
    1.60660297316774, 2.40965937926238, 2.39721551467612, 2.45333816867543,
    2.17369588423621, 1.83988225331191, 2.62537066751224, 1.40236145181362,
    2.02069874809587, 2.25007385714502, 1.88697660890813, 2.71687671988929,
    2.72541641338052, 3.84984414187981, 0.990461286855624, 0.992006460987693,
    0.823798773849397, 0.670068730341277, 0.515036259090417, 0.65575291763732,
    0.799234706529325, 0.763770392688469, 0.85039662021339, 0.856659788271371,
    0.871545121350289, 0.893180671332443, 0.77003356074645, 0.866744655894669,
    0.846329062125548, 0.770277863202376, 0.764826961908685, 0.666082087716545,
    0.991844630061474, 0.991600327605548, 0.991763714598365, 0.991925545524583,
    0.989404717635408, 0.989649020091335, 0.950036232699527, 0.469650464619079,
    0.678283205913471, 0.605160524328217, 0.715373609349505, 0.740101063662394,
    0.590763796161151, 0.531141548382368, 0.730665075810569, 0.75238309732243,
    0.755555917116275, 0.762469521012327, 0.723670556451407, -0.521707116597141,
    -0.422635468419367, 0.619882470414319, 0.716676037092246, 0.736114421037662,
    0.742866194007495, 0.739612458750542, 0.737905453692255, 0.725460033039402,
    0.741401935338536, 0.748804144146441, 0.743518185912165, 0.749292749058294,
    0.712689394467514, 0.755229143130641, 0.764339913063431, 0.673322465598231,
    0.792889066941637, 0.821928381798293, 0.807855315483664, 0.810622001895364,
    0.643957932822539, 0.755719304109091, 0.653963441049326, 0.776702862186576,
    0.746200844727558, 0.780770420274418, 0.835836505053507, 0.858285877866549,
    0.85112797151457, 0.852674701713237, 0.859018785234327, 0.79890948861029,
    0.855358916595229, 0.822742204629181, 0.86446813046142, 0.717651690849353,
    0.812086260564323, 0.799396537455544, 0.631677455229101, 0.759704390667226,
    0.76149386725522, 0.537404716440349, 0.490960796682197, -0.00463552239620513,
    0.115663986314979, 0.256053870876249, 0.595724536476392, 0.59792170251313,
    0.719196864981421, 0.7833737196933, 0.39774618318026, -0.66494460303322,
    -0.814526172990389, -0.590763796161151, -0.605404826784144,
    -0.447282007269146, -0.424832634456105, -0.459970174311327,
    -0.355287349980242, -0.416291384898277, -0.328445201160323,
    -0.361875735957258, -0.28696357778401, -0.374321156610111,
    -0.437033752653031, -0.50633317860637, -0.570754335774175,
    -0.653474836137473, -0.417755643567236, -0.5866153226102, -0.379689586374095,
    -0.392866358328128, -0.427028244426244, -0.578562677964225,
    -0.280862240652247, -0.780364286892273, -0.775240159584215,
    -0.469487077626262, -0.465582906531237, -0.49144940159405,
    -0.426783941970318, -0.443133533718196, -0.459970174311327,
    -0.532686722514436, -0.451430480820097, -0.428492503095203,
    -0.317952644088281, -0.495109270233148, -0.616384432701439,
    -0.50975030085614, -0.407507388951121, -0.358460169774086,
    -0.388473582321251, -0.421415512206334, -0.39555057321012,
    -0.615653081400259, -0.566850164679151, -0.418731297324342,
    -0.479248283397123, -0.523414121655428, -0.632734024449317,
    -0.471439941207074, -0.491205099138123, -0.562945993584126,
    -0.398479090548038, -0.440205016380277, -0.584175410184135,
    -0.473391248721287, -0.563188739973454, -0.43508088907222,
    -0.511701608370354, -0.76572325626928, -0.815746129203422,
    -0.519754253016329, -0.477295419816311, -0.467291467656123,
    -0.448013358570327, -0.542692230741224, -0.683001199839384,
    -0.455578954371048, -0.673728598980376, -0.600523445865414,
    -0.468511423869155, -0.488032279344279, -0.58807958127916,
    -0.375296810367218, -0.558797520033175, -0.44923331478336, -0.52951545878719,
    -0.444842094843081, -0.641518020396473, -0.448746265938106,
    -0.699595094043188, -0.483396756948074, -0.524147029023207,
    -0.58246684905925, -0.399454744305145, -0.564652998642413,
    -0.468022818957303, -0.644935142646243, -0.518778599259223,
    -0.50072044638646, -0.542203625829371, -0.611748910305234,
    -0.800860796124503, -0.441913577505163, -0.461434432980286,
    -0.418731297324342, -0.384813713682153, -1.16615365433153, -1.99434364812121,
    -1.9143073626399, -2.20785932640093, -1.92309135858706, -2.08584970049826,
    -2.03314261268213, -0.54928061671824, -4.02846347062704, -1.79888767878085,
    -1.01022644478667, -2.123917313758, -2.91428555281046, -1.62197535509191,
    -2.29131117806541, -1.00925079102957, -2.58339888315747, -1.72202265702679,
    -2.01606166963307, -1.92797118343919, -2.19541390574808, -2.79959722025259,
    -3.31252034090257, -0.748398010764296, -1.38991603116077, -1.37063947814157,
    -2.07145297233119, -1.63710499062675, -1.97506709510201, -2.01801453321388,
    -2.77910071102036, -1.882096784056, -3.42915998097466, -1.50362715389494,
    -1.4858148595447, -2.91794542144956, -2.26495763415734, -1.61245845177697,
    -2.65025995275135, -1.08538446148245, -1.34819010532853, -1.80962453830882,
    -1.78888217055406, -1.29719002257068, -1.92845978835104, -2.07950561697717,
    -2.58535174673829, -2.73493331669545, -1.16054092211162, -2.07120866987527,
    -2.21249484879713, -1.06073792263267, -2.59120722534752, -3.14585627182975,
    -1.3396488557707, -3.59313827909889, -3.80811665998098, -4.69535938117447,
    48.6406788834742, 48.645966397775, 48.7781511431637, 42.7753488166026,
    39.1739571602666, 41.2396355694517, 43.0631231050843, 43.4198786039,
    44.4458812298931, 44.3783697244612, 44.1473700818835, 44.6612991995892,
    44.5912785367676, 44.8559475703651, 44.8875497269091, 44.9315132765096,
    44.0463868058897, 46.7308436562377, 48.6369007537736, 48.6425912893234,
    48.6373162235553, 48.6375566358448, 48.6580484768772, 48.6258830242247,
    49.1398985036266, 35.7881507507859, 35.9714810711594, 37.358634306075,
    38.0181739117791, 37.9959221594241, 36.3282785837733, 38.1011745041301,
    37.4605527780943, 37.7191656004913, 37.7230759958528, 37.6982645139436,
    38.1330606428283, 10.1452601240432, 21.874757098004, 36.5271952452281,
    36.5370653756609, 36.9320954445939, 36.8900108454107, 36.8027248456482,
    36.7803120404002, 36.7779118076724, 37.1248267413207, 37.268514709037,
    37.3083266729533, 37.4749938541594, 37.4391553063016, 37.7397056795883,
    38.203573022695, 38.5135687202331, 40.7760374258665, 40.6764849531099,
    40.7072701746898, 40.6826220797734, 39.6924447750719, 39.1258062354502,
    38.8272025014779, 39.9459754839591, 40.3183857907866, 41.2161591926838,
    42.2207301134191, 42.6059133929183, 42.6415053042202, 42.6814051858994,
    42.7971882113092, 42.5700834034269, 42.4885120582464, 42.7285765667923,
    42.765499693069, 42.7787254811209, 42.6506129620198, 42.721053762823,
    38.8437963956817, 38.1541702423002, 38.1819452530458, 35.4977661605857,
    33.7359346517757, 39.1241326858237, 26.9684267365006, 32.2401096158044,
    32.5978353221441, 35.1565977808847, 38.5611890263395, 39.1479677359425,
    37.9059566129793, 1.488492072127, 2.66025923671174, 2.74432495664923,
    2.42356516841769, 3.19075501745583, 2.96247849142505, 3.06325947876101,
    1.88685523571347, 3.89583991445767, 2.76140434363169, 3.72685652812014,
    2.50445962462794, 3.24065885129521, 2.30802411136398, 3.85948008428951,
    4.34458618039577, 2.68490888769472, 2.95686653723844, 3.06301906647158,
    3.9724614117891, 3.57374113091833, 2.20541941397486, 2.35475512540951,
    4.01382088393745, 2.49152793316313, 2.43454866450149, 3.73649791676294,
    2.96418705254994, 3.36217831621942, 2.21835421757287, 3.12353294041116,
    3.560563580931, 2.60145470193034, 3.59667677454333, 2.94344001659539,
    3.60338653371501, 2.4457787971411, 1.94212983341672, 3.63498791222569,
    3.15061705758711, 3.5067524638634, 3.15683665578034, 3.38365203527536,
    4.00577057339137, 3.70404070162321, 3.16525264197701, 3.17783110632401,
    4.14912943302211, 3.40013000251759, 3.13255734864775, 2.58315535873485,
    3.12792571641804, 3.65145498700174, 2.05388186830369, 4.00320928777064,
    2.50323500021511, 3.30508156452961, 3.98942253771009, 3.54884951157933,
    3.11608560567191, 2.25495679413035, 1.72800262096389, 3.5817984437641,
    4.03773140328625, 2.49042701604486, 2.41710126776854, 4.54674801874132,
    2.31741886345095, 3.57874466306502, 2.95625889323182, 3.11437860061362,
    3.23748680953466, 3.49906938503451, 2.81765770722525, 2.60242257535446,
    3.83130127426169, 3.85069686637566, 3.35083381268484, 3.47637337566565,
    3.87766505659003, 3.37096075609998, 2.60352893870582, 2.76531162685991,
    3.24663025686601, 4.86909888506957, 2.647091023124, 3.8299568327208,
    3.95635223233065, 4.33165682303086, 3.45929710081639, 3.48211292731332,
    3.4791812978422, 3.97746727803569, 2.81326415318507, 3.79347874149118,
    4.4647588696792, 4.1213038501121, 3.26249746796843, 3.2971471809451,
    1.8053562476298, 7.90258380750625, 6.07222156964247, 6.61247077582443,
    7.02535359683923, 3.19587914476389, 16.4388788492855, 5.2702629685329,
    8.44857178730258, 7.34806936078477, 8.19662437617913, 4.16401398806773,
    7.0613407490571, 6.41164637672009, 9.50760737307675, 7.85853856440931,
    8.2667819728614, 2.88841983578095, 1.85659518661048, 10.4465791942294,
    9.15792341891683, 3.37315714410342, 9.01358968356228, 9.56275281725235,
    8.95930396618222, 5.62762844545511, 7.15419591117484, 3.30775644301201,
    4.94230715031773, 5.53794921526441, 9.51455910060446, 5.91970992628079,
    9.50430773385515, 9.76260311866619, 7.40637828835463, 9.5479888573681,
    -3.93402812287877, 1.25741773835256, 8.04752518686742, 6.30049576157335,
    1.94114873342652, 11.7917178529587, 6.15470321378294, 14.1407365388538,
    9.18086061860842, -0.770352554399093, -1.39308807292132, 18.7566696126534,
    -1.51033380093342, 11.1685015098576, 6.78035946931009, 7.91990827497793,
    8.94258792075045, 17.055141130759, 20.6176653556038, 2.75897298957191,
    24.1274677099998, 16.0, 17.0, 16.0, 13.0, 16.0, 17.0, 16.0, 17.0, 16.0, 17.0,
    16.0, 17.0, 18.0, 16.0, 19.0, 19.0, 16.0, 15.0, 17.0, 15.0, 15.0, 16.0, 15.0,
    14.0, 14.0, 16.0, 16.0, 15.0, 14.0, 18.0, 18.0, 17.0, 21.0, 21.0, 21.0, 20.0,
    19.0, 7.0, 10.0, 16.0, 19.0, 21.0, 18.0, 19.0, 22.0, 19.0, 18.0, 19.0, 20.0,
    19.0, 16.0, 17.0, 20.0, 17.0, 16.0, 18.0, 21.0, 19.0, 17.0, 17.0, 17.0, 16.0,
    18.0, 13.0, 18.0, 20.0, 17.0, 21.0, 18.0, 20.0, 20.0, 20.0, 15.0, 19.0, 18.0,
    18.0, 15.0, 18.0, 18.0, 14.0, 14.0, 18.0, 12.0, 14.0, 15.0, 13.0, 18.0, 17.0,
    18.0, 14.0, 13.0, 12.0, 12.0, 12.0, 10.0, 11.0, 11.0, 13.0, 11.0, 14.0, 11.0,
    11.0, 12.0, 11.0, 15.0, 10.0, 11.0, 10.0, 10.0, 10.0, 10.0, 9.0, 13.0, 13.0,
    9.0, 16.0, 10.0, 10.0, 11.0, 13.0, 14.0, 11.0, 12.0, 15.0, 13.0, 12.0, 11.0,
    15.0, 14.0, 13.0, 12.0, 14.0, 14.0, 11.0, 13.0, 13.0, 12.0, 14.0, 13.0, 12.0,
    13.0, 14.0, 10.0, 11.0, 12.0, 13.0, 14.0, 14.0, 12.0, 12.0, 12.0, 15.0, 13.0,
    12.0, 14.0, 13.0, 13.0, 14.0, 10.0, 12.0, 12.0, 11.0, 13.0, 14.0, 11.0, 11.0,
    12.0, 12.0, 14.0, 9.0, 14.0, 15.0, 13.0, 14.0, 12.0, 13.0, 12.0, 14.0, 12.0,
    16.0, 14.0, 10.0, 13.0, 11.0, 13.0, 16.0, 11.0, 13.0, 15.0, 18.0, 15.0, 18.0,
    17.0, 17.0, 15.0, 16.0, 15.0, 15.0, 14.0, 17.0, 14.0, 19.0, 16.0, 15.0, 17.0,
    18.0, 18.0, 17.0, 15.0, 15.0, 18.0, 15.0, 14.0, 16.0, 16.0, 18.0, 16.0, 15.0,
    15.0, 16.0, 17.0, 15.0, 16.0, 17.0, 12.0, 15.0, 15.0, 16.0, 14.0, 17.0, 15.0,
    14.0, 16.0, 17.0, 15.0, 15.0, 12.0, 15.0, 16.0, 15.0, 14.0, 16.0, 14.0, 15.0,
    13.0, 43.0, 46.0, 44.0, 45.0, 44.0, 44.0, 44.0, 45.0, 46.0, 43.0, 44.0, 46.0,
    43.0, 45.0, 46.0, 46.0, 44.0, 41.0, 47.0, 43.0, 42.0, 46.0, 46.0, 43.0, 45.0,
    45.0, 45.0, 47.0, 46.0, 46.0, 45.0, 43.0, 45.0, 47.0, 45.0, 43.0, 44.0, 18.0,
    30.0, 44.0, 44.0, 47.0, 46.0, 45.0, 46.0, 44.0, 44.0, 46.0, 44.0, 44.0, 44.0,
    45.0, 46.0, 44.0, 45.0, 43.0, 47.0, 47.0, 42.0, 46.0, 47.0, 44.0, 43.0, 44.0,
    46.0, 45.0, 44.0, 46.0, 44.0, 45.0, 47.0, 46.0, 45.0, 47.0, 43.0, 44.0, 46.0,
    45.0, 45.0, 46.0, 40.0, 42.0, 42.0, 46.0, 45.0, 42.0, 46.0, 45.0, 44.0, 44.0,
    42.0, 46.0, 42.0, 42.0, 45.0, 45.0, 44.0, 46.0, 39.0, 45.0, 40.0, 44.0, 44.0,
    43.0, 42.0, 38.0, 44.0, 40.0, 42.0, 39.0, 39.0, 39.0, 44.0, 43.0, 36.0, 43.0,
    39.0, 42.0, 43.0, 46.0, 45.0, 44.0, 43.0, 45.0, 46.0, 43.0, 41.0, 46.0, 47.0,
    45.0, 41.0, 45.0, 44.0, 46.0, 41.0, 45.0, 46.0, 43.0, 43.0, 41.0, 45.0, 43.0,
    43.0, 46.0, 43.0, 46.0, 42.0, 46.0, 47.0, 44.0, 41.0, 47.0, 45.0, 42.0, 45.0,
    46.0, 45.0, 46.0, 39.0, 41.0, 44.0, 43.0, 44.0, 45.0, 44.0, 44.0, 46.0, 38.0,
    44.0, 39.0, 45.0, 46.0, 45.0, 41.0, 44.0, 47.0, 44.0, 47.0, 40.0, 46.0, 45.0,
    43.0, 41.0, 40.0, 42.0, 46.0, 47.0, 40.0, 46.0, 46.0, 44.0, 43.0, 45.0, 45.0,
    47.0, 47.0, 43.0, 46.0, 45.0, 43.0, 46.0, 47.0, 45.0, 43.0, 42.0, 46.0, 45.0,
    45.0, 40.0, 40.0, 46.0, 42.0, 42.0, 44.0, 46.0, 45.0, 43.0, 44.0, 46.0, 46.0,
    47.0, 43.0, 45.0, 46.0, 45.0, 44.0, 45.0, 43.0, 43.0, 46.0, 45.0, 41.0, 43.0,
    46.0, 44.0, 47.0, 40.0, 41.0, 45.0, 43.0, 47.0, 44.0, 44.0, 42.0, 41.0,
    2.6875, 2.70588235294118, 2.75, 3.46153846153846, 2.75, 2.58823529411765,
    2.75, 2.64705882352941, 2.875, 2.52941176470588, 2.75, 2.70588235294118,
    2.38888888888889, 2.8125, 2.42105263157895, 2.42105263157895, 2.75,
    2.73333333333333, 2.76470588235294, 2.86666666666667, 2.8, 2.875,
    3.06666666666667, 3.07142857142857, 3.21428571428571, 2.8125, 2.8125,
    3.13333333333333, 3.28571428571429, 2.55555555555556, 2.5, 2.52941176470588,
    2.14285714285714, 2.23809523809524, 2.14285714285714, 2.15, 2.31578947368421,
    2.57142857142857, 3.0, 2.75, 2.31578947368421, 2.23809523809524,
    2.55555555555556, 2.36842105263158, 2.09090909090909, 2.31578947368421,
    2.44444444444444, 2.42105263157895, 2.2, 2.31578947368421, 2.75,
    2.64705882352941, 2.3, 2.58823529411765, 2.8125, 2.38888888888889,
    2.23809523809524, 2.47368421052632, 2.47058823529412, 2.70588235294118,
    2.76470588235294, 2.75, 2.38888888888889, 3.38461538461538, 2.55555555555556,
    2.25, 2.58823529411765, 2.19047619047619, 2.44444444444444, 2.25, 2.35, 2.3,
    3.0, 2.47368421052632, 2.38888888888889, 2.44444444444444, 3.06666666666667,
    2.5, 2.5, 3.28571428571429, 2.85714285714286, 2.33333333333333, 3.5,
    3.28571428571429, 3.0, 3.23076923076923, 2.55555555555556, 2.64705882352941,
    2.44444444444444, 3.14285714285714, 3.23076923076923, 3.83333333333333, 3.5,
    3.5, 4.5, 4.09090909090909, 4.0, 3.53846153846154, 3.54545454545455,
    3.21428571428571, 3.63636363636364, 4.0, 3.66666666666667, 3.90909090909091,
    2.8, 3.8, 4.0, 4.0, 4.2, 3.9, 3.9, 4.33333333333333, 3.38461538461538,
    3.30769230769231, 4.0, 2.6875, 3.9, 4.2, 3.90909090909091, 3.53846153846154,
    3.21428571428571, 4.0, 3.58333333333333, 3.0, 3.53846153846154,
    3.58333333333333, 3.72727272727273, 3.06666666666667, 3.35714285714286,
    3.46153846153846, 3.41666666666667, 3.21428571428571, 3.14285714285714,
    4.18181818181818, 3.15384615384615, 3.46153846153846, 3.83333333333333,
    3.07142857142857, 3.30769230769231, 3.41666666666667, 3.46153846153846,
    3.07142857142857, 4.3, 4.18181818181818, 3.58333333333333, 3.53846153846154,
    3.0, 3.28571428571429, 3.91666666666667, 3.66666666666667, 3.41666666666667,
    3.13333333333333, 3.46153846153846, 3.5, 3.21428571428571, 3.53846153846154,
    3.46153846153846, 3.28571428571429, 3.9, 3.41666666666667, 3.66666666666667,
    3.90909090909091, 3.38461538461538, 3.21428571428571, 4.0, 4.0,
    3.83333333333333, 3.16666666666667, 3.14285714285714, 4.33333333333333,
    3.21428571428571, 3.06666666666667, 3.46153846153846, 2.92857142857143,
    3.66666666666667, 3.61538461538462, 3.66666666666667, 3.35714285714286,
    3.33333333333333, 2.875, 3.21428571428571, 4.3, 3.15384615384615,
    3.63636363636364, 3.23076923076923, 2.875, 4.27272727272727,
    3.07692307692308, 3.06666666666667, 2.55555555555556, 2.93333333333333,
    2.38888888888889, 2.64705882352941, 2.64705882352941, 3.13333333333333,
    2.9375, 2.86666666666667, 3.06666666666667, 3.21428571428571,
    2.52941176470588, 3.28571428571429, 2.47368421052632, 2.8125,
    2.86666666666667, 2.47058823529412, 2.55555555555556, 2.5, 2.64705882352941,
    2.66666666666667, 2.66666666666667, 2.55555555555556, 2.8, 3.0, 2.75, 2.875,
    2.5, 2.6875, 2.93333333333333, 3.06666666666667, 2.875, 2.76470588235294,
    2.86666666666667, 2.8125, 2.70588235294118, 3.75, 2.93333333333333, 3.0,
    2.6875, 3.07142857142857, 2.70588235294118, 3.0, 2.92857142857143, 2.6875,
    2.70588235294118, 2.93333333333333, 3.13333333333333, 3.33333333333333,
    2.73333333333333, 2.8125, 2.86666666666667, 3.35714285714286, 2.75,
    3.14285714285714, 2.8, 3.15384615384615, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 2.0, 3.0, 2.0, 2.0, 3.0, 4.0, 2.0, 5.0, 3.0, 4.0, 5.0, 3.0, 3.0,
    3.0, 6.0, 4.0, 2.0, 3.0, 4.0, 4.0, 1.0, 2.0, 4.0, 4.0, 4.0, 2.0, 3.0, 4.0,
    2.0, 5.0, 3.0, 2.0, 5.0, 4.0, 3.0, 1.0, 2.0, 3.0, 2.0, 2.0, 3.0, 2.0, 5.0,
    8.0, 5.0, 4.0, 7.0, 4.0, 7.0, 3.0, 3.0, 5.0, 8.0, 10.0, 6.0, 14.0, 50.0,
    50.0, 49.0, 49.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 49.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 19.0, 35.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 47.0, 42.0, 49.0, 50.0, 50.0, 50.0, 50.0, 49.0, 1.0, 5.0, 1.0,
    2.0, 2.0, 2.0, 1.0, 0.0, 3.0, 1.0, 2.0, 0.0, 2.0, 2.0, 3.0, 1.0, 0.0, 0.0,
    1.0, 1.0, 0.0, 2.0, 0.0, 1.0, 4.0, 1.0, 5.0, 4.0, 3.0, 1.0, 2.0, 3.0, 2.0,
    6.0, 2.0, 3.0, 4.0, 2.0, 5.0, 1.0, 4.0, 3.0, 4.0, 0.0, 6.0, 4.0, 5.0, 2.0,
    1.0, 2.0, 3.0, 2.0, 5.0, 4.0, 5.0, 5.0, 5.0, 4.0, 3.0, 3.0, 4.0, 2.0, 5.0,
    3.0, 4.0, 5.0, 4.0, 5.0, 3.0, 4.0, 3.0, 5.0, 4.0, 4.0, 2.0, 4.0, 3.0, 4.0,
    6.0, 6.0, 6.0, 3.0, 5.0, 3.0, 6.0, 4.0, 4.0, 3.0, 5.0, 8.0, 6.0, 4.0, 7.0,
    5.0, 6.0, 5.0, 3.0, 5.0, 5.0, 10.0, 13.0, 10.0, 13.0, 12.0, 11.0, 10.0, 14.0,
    13.0, 11.0, 13.0, 11.0, 11.0, 13.0, 15.0, 13.0, 13.0, 12.0, 15.0, 10.0, 9.0,
    13.0, 12.0, 12.0, 10.0, 11.0, 12.0, 10.0, 12.0, 13.0, 8.0, 11.0, 14.0, 12.0,
    11.0, 13.0, 12.0, 14.0, 12.0, 11.0, 7.0, 18.0, 10.0, 17.0, 10.0, 9.0, 12.0,
    12.0, 10.0, 11.0, 16.0, 15.0, 11.0, 11.0, 13.0, 9.0, 8.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 31.0, 15.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 3.0, 8.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 49.0, 45.0, 49.0, 48.0,
    48.0, 48.0, 49.0, 50.0, 47.0, 49.0, 48.0, 50.0, 48.0, 48.0, 47.0, 49.0, 50.0,
    50.0, 49.0, 49.0, 50.0, 48.0, 50.0, 49.0, 46.0, 49.0, 45.0, 46.0, 47.0, 49.0,
    48.0, 47.0, 48.0, 44.0, 48.0, 47.0, 46.0, 48.0, 45.0, 49.0, 46.0, 47.0, 46.0,
    50.0, 44.0, 46.0, 45.0, 48.0, 49.0, 48.0, 47.0, 48.0, 45.0, 46.0, 45.0, 45.0,
    45.0, 46.0, 47.0, 47.0, 46.0, 48.0, 45.0, 47.0, 46.0, 45.0, 46.0, 45.0, 47.0,
    46.0, 47.0, 45.0, 46.0, 46.0, 48.0, 46.0, 47.0, 46.0, 44.0, 44.0, 44.0, 47.0,
    45.0, 47.0, 44.0, 46.0, 46.0, 47.0, 45.0, 42.0, 44.0, 46.0, 43.0, 45.0, 44.0,
    45.0, 47.0, 45.0, 45.0, 39.0, 35.0, 37.0, 35.0, 36.0, 36.0, 36.0, 34.0, 32.0,
    36.0, 33.0, 34.0, 36.0, 34.0, 32.0, 31.0, 33.0, 36.0, 32.0, 36.0, 37.0, 36.0,
    36.0, 34.0, 36.0, 35.0, 36.0, 37.0, 34.0, 35.0, 37.0, 36.0, 34.0, 33.0, 35.0,
    34.0, 37.0, 34.0, 35.0, 37.0, 41.0, 29.0, 38.0, 28.0, 32.0, 36.0, 34.0, 31.0,
    36.0, 32.0, 31.0, 32.0, 34.0, 31.0, 27.0, 35.0, 28.0, 0.000830916698299489,
    1.40377232125122E-5, 0.104338201240155, 0.00579690655918919,
    0.105059633932124, 0.0755817596350952, 0.0196482333446578,
    0.0468413079720106, 0.0120010223036061, 0.0120184158873187,
    0.00531599039814257, 0.000274840110535801, 0.00279081527320257,
    0.00900680888547489, 0.0025114597049154, 0.099681756067018,
    0.0144947802626274, 0.319121235325986, 1.93550170063395E-5,
    0.000698263994519222, 0.000326892650873489, 0.000259299574123406,
    0.0010432323338293, 0.00056878800975646, 0.0144002884322239,
    0.0718891402764003, 0.00243063764088936, 0.0195837357921315,
    0.0153024600303506, 0.000876930670968168, 0.0219711828817337,
    0.117302194070274, 0.00248090194966244, 0.00471727432672031,
    0.0037336538163586, 0.00106695040438709, 0.0110487473545406,
    0.0152261423962767, 0.166947998845689, 0.0366789237157232,
    0.00718185812593386, 4.78906949585856E-6, 0.000790621441139283,
    0.0025477836420893, 0.00141084951369519, 0.000256035975985336,
    0.00419932766310704, 0.000692447550276532, 0.00305511153037674,
    0.000548847307186851, 0.00489065339575601, 0.00139090458887081,
    0.000444125577126475, 0.00487631770787622, 0.0163801567504664,
    0.000472326701945858, 0.00235564478046624, 0.00115502792942912,
    0.0232245444530494, 0.0440001953823682, 0.0179770170007815,
    0.00283293436718335, 0.023427044726231, 0.0447836650359056,
    0.00198569122297254, 0.003616503237585, 0.00599881270388258,
    0.00387076545474607, 0.00231439111415233, 0.0138856926481794,
    0.000644860466021857, 8.5172333680994E-5, 7.62482575827894E-5,
    0.0232728892196477, 0.0524238898529948, 0.028477814389936,
    0.0953231835433056, 0.00823625154400554, 0.0118623708929126,
    0.105403085992518, 0.0201570534775017, 0.160089216050703, 0.0225958927670071,
    0.133972270268411, 0.0030379659783526, 0.00548705698230712,
    0.00535164686424172, 0.00712826491802909, 0.200197522120997,
    0.176871404386832, 0.100024165368051, 0.222854715125583, 0.204416703263549,
    0.101596677188847, 0.330237843816716, 0.0499396049141244, 0.296335299073888,
    0.174453534304139, 0.257735596335941, 0.121514158363649, 0.24121360973132,
    0.205307080451163, 0.117646974755339, 0.109376752467361, 0.367784773062382,
    0.247935379323541, 0.0396185755485088, 0.283625040508275, 0.227524981219249,
    0.100132235334203, 0.0484700089511709, 0.100174407524714, 0.137863474877467,
    0.393783154951385, 0.237397843011627, 0.279427159458114, 0.231829855984741,
    0.596960105023992, 0.112226516191482, 0.254955152480376, 0.323279338653107,
    0.306321292955311, 0.323391192609369, 0.351137987874426, 0.195400219255149,
    0.0281108488921304, 0.0727033136190497, 0.358127867146119,
    0.0538413989932989, 0.293830635669339, 0.150921979064238, 0.0497641351113073,
    0.17548116369631, 0.0265529554397619, 0.0249134866579908, 0.50821817568492,
    0.42789325632794, 0.168704837847261, 0.265516119963293, 0.139836305633881,
    0.0809608054215804, 0.0788146635038614, 0.0626640715629515,
    0.0745881176505774, 0.292521506489047, 0.0365180052129517,
    0.0797760427385773, 0.0355453742009262, 0.0804414844234297,
    0.253190397828498, 0.0253384687179085, 0.293495385028231, 0.141938929976558,
    0.442753582077018, 0.150835171735331, 0.0360448959528167, 0.249834046237669,
    0.159112631816974, 0.230693896009608, 0.151128816328182, 0.160642664691591,
    0.0319258369124902, 0.452927768800945, 0.0110648524749811, 0.234424254899122,
    0.2421542576336, 0.105155828326595, 0.223944443989654, 0.250094496580429,
    0.273602066313867, 0.313960669084357, 0.0106942191777222, 0.125681946487929,
    0.324544114022737, 0.199109106424872, 0.0817971467599115, 0.0201627501279178,
    0.26832772020141, 0.311138432329372, 0.3174527606891, 0.120717246960652,
    0.139474740750324, 0.0824544220391772, 0.358742147400049, 0.335208435140226,
    0.183241957424814, 0.221247570979969, 0.274407262019165, 1.5397795337181,
    1.66182078218415, 2.16269655473493, 2.18111011672142, 2.06664726166585,
    2.07578468914207, 2.19856890378892, 1.53711860843728, 1.38837303485216,
    2.45555084270727, 2.12604530900617, 1.03827404125579, 1.83085116937932,
    0.990966317094369, 1.36580561937285, 1.89910332599439, 2.21812651000093,
    1.28559902156472, 0.982176626861627, 2.61107794759839, 2.00634013533084,
    3.11116659165789, 2.1940853248263, 2.86063689565166, 2.45424564941686,
    1.77352780338008, 2.64278973838139, 1.87652334183264, 1.86193542022079,
    2.50756182643117, 2.60988158428341, 2.60963208178206, 2.58858585994506,
    2.2410798742958, 2.01169677913817, 2.57747600221102, 0.630718424494917,
    2.52144906777289, 2.44958459308265, 2.77017409656782, 1.93527420851852,
    2.01662723874892, 2.11015946886152, 2.39575712176997, 2.02536296833141,
    1.93314003720381, 1.60033273336522, 2.54622745068797, 2.24671129540925,
    2.02115895301093, 2.13199934571098, 2.8451414226137, 2.55712128901748,
    2.4259987160104, 0.248944278367062, 0.0330625896000125, 2.85348996623178,
    0.000830916698299489, 2.80754464250244E-6, -0.104338201240155,
    0.0019323021863964, 0.00750425956658032, -0.00419898664639418,
    -0.0196482333446578, -0.0468413079720106, -0.00600051115180306,
    0.00042922913883281, -0.00531599039814257, -0.000274840110535801,
    0.00279081527320257, 0.00900680888547489, -0.0025114597049154,
    -0.002431262343098, 0.0144947802626274, 0.319121235325986,
    -1.07527872257441E-6, 0.000139652798903844, 0.000163446325436745,
    -1.85213981516719E-5, -0.0010432323338293, 0.00056878800975646,
    -0.00480009614407462, -0.0179722850691001, 0.000347233948698479,
    -0.0195837357921315, -0.000493627742914536, -3.37281027295449E-5,
    -0.0219711828817337, -0.117302194070274, 0.000496180389932488,
    0.00471727432672031, 0.0037336538163586, -0.00106695040438709,
    -0.0110487473545406, 0.00190326779953459, 0.166947998845689,
    -0.00152828848815513, -0.00102597973227627, 4.78906949585856E-6,
    -0.000790621441139283, 0.0025477836420893, 0.000128259046699563,
    0.000256035975985336, -0.000127252353427486, 0.000692447550276532,
    0.00305511153037674, -0.000548847307186851, 0.00489065339575601,
    -0.000126445871715528, -3.70104647605396E-5, -0.00487631770787622,
    0.00819007837523322, 6.74752431351225E-5, -0.00235564478046624,
    -0.000231005585885824, 0.00580613611326236, 0.0440001953823682,
    -0.00163427427279832, -0.00283293436718335, 0.023427044726231,
    -0.0447836650359056, 0.000248211402871568, -0.003616503237585,
    0.00599881270388258, 0.00387076545474607, -0.000330627302021761,
    -0.0138856926481794, 0.000644860466021857, 8.5172333680994E-5,
    -6.35402146523245E-6, 0.0232728892196477, -0.0524238898529948,
    -0.000889931699685501, -0.0953231835433056, -0.000588303681714681,
    -0.0118623708929126, 0.00234229079983374, 0.00335950891291696,
    -0.160089216050703, -0.0225958927670071, 0.0148858078076012,
    0.000759491494588149, 0.00548705698230712, 0.00267582343212086,
    -0.00712826491802909, -0.200197522120997, -0.176871404386832,
    -0.0111137961520056, -0.0117291955329254, -0.204416703263549,
    -0.101596677188847, -0.0235884174154797, -0.00226998204155111,
    0.00779829734404969, 0.0158594122094672, -0.011205895492867,
    -0.00391981156011771, -0.24121360973132, -0.0205307080451163,
    0.0117646974755339, -0.00437507009869445, -0.0334349793693075,
    0.00918279182679781, 0.000900422171557018, -0.011345001620331,
    -0.227524981219249, 0.0055629019630113, -0.0484700089511709,
    0.00556524486248412, -0.0137863474877467, -0.393783154951385,
    0.00719387403065537, -0.279427159458114, 0.0136370503520436,
    0.596960105023992, -0.00510120528143102, -0.254955152480376,
    0.00873727942305695, 0.0105628032053555, -0.0153995806004461,
    0.00949021588849799, -0.0977001096275746, 0.00117128537050544,
    0.00454395710119061, -0.358127867146119, -0.0538413989932989,
    -0.293830635669339, 0.00486845093755605, 0.00261916500585828,
    -0.17548116369631, -0.0265529554397619, -0.0249134866579908,
    -0.50821817568492, -0.42789325632794, 0.00992381399101537,
    0.0265516119963293, 0.0139836305633881, 0.00449782252342113,
    -0.00414814018441376, 0.00174066865452643, -0.0745881176505774,
    -0.292521506489047, -0.0365180052129517, -0.00797760427385773,
    -0.0355453742009262, 0.00217409417360621, -0.253190397828498,
    0.000745249079938484, -0.293495385028231, 0.00788549610980877,
    0.0442753582077018, 0.0100556781156887, 0.00400498843920185,
    0.0624585115594172, -0.0079556315908487, 0.00961224566706701,
    -0.00755644081640911, -0.160642664691591, -0.0319258369124902,
    -0.452927768800945, -0.000526897736903863, -0.234424254899122,
    -0.0269060286259555, -0.105155828326595, -0.223944443989654,
    0.0625236241451073, -0.00855006457230833, 0.018468274652021,
    0.000629071716336599, -0.00661483928883838, -0.0170812691590914,
    -0.199109106424872, -0.0102246433449889, 0.000438320654954735,
    -0.26832772020141, -0.311138432329372, -0.3174527606891, 0.0134130274400724,
    -0.139474740750324, 0.00916160244879747, -0.358742147400049,
    -0.335208435140226, -0.00916209787124069, -0.221247570979969,
    0.0274407262019165, -0.10265196891454, -0.11078805214561,
    -0.0800998723975899, -2.18111011672142, -2.06664726166585,
    -0.0518946172285518, -2.19856890378892, -0.0569303188310104,
    -0.0694186517426079, -2.45555084270727, -2.12604530900617, -1.03827404125579,
    -1.83085116937932, 0.090087847008579, 0.0426814256054016, -1.89910332599439,
    -2.21812651000093, -0.0476147785764713, 0.0892887842601479,
    -0.186505567685599, 0.0542254090629957, -3.11116659165789, -2.1940853248263,
    -2.86063689565166, -2.45424564941686, -1.77352780338008, -2.64278973838139,
    -0.06950086451232, 0.0517204283394665, -2.50756182643117, -2.60988158428341,
    -2.60963208178206, -2.58858585994506, -2.2410798742958, -2.01169677913817,
    -2.57747600221102, 0.0630718424494917, 0.252144906777289, -2.44958459308265,
    -1.38508704828391, -1.93527420851852, -2.01662723874892, -2.11015946886152,
    -0.171125508697855, 0.0880592594926701, -1.93314003720381,
    -0.228618961909317, -0.19586365005292, 0.0680821604669469, -2.02115895301093,
    -2.13199934571098, 0.129324610118805, -2.55712128901748, -2.4259987160104,
    -0.0207453565305885, -0.00110208632000042, -2.85348996623178,
    -0.349774091516507, -0.343993243480765, 0.00490873862807961,
    -0.0859928843431683, -0.0737418280459186, 0.339239505623672,
    -0.125300814222444, -0.0825422614103153, -0.177933449071352,
    -0.201211035026315, -0.178081792781984, -0.19610804063227,
    -0.184986131273415, -0.185292765247501, -0.262716757517453,
    -0.186201385625803, -0.18566137680428, -0.534858228877706,
    -0.358953582472775, -0.353366111154788, -0.358052875817821,
    -0.358313986001689, -0.337930711338224, -0.369648358139721,
    0.297879023049356, -0.578598385383647, -0.108413551032059,
    -0.123602507364043, -0.163487846100769, -0.165150809030357,
    -0.120446415023898, -0.0706186298426416, -0.168999857309767,
    -0.178287001438906, -0.174263843989894, -0.171591675339964,
    -0.131056946885395, 0.356243975235761, 0.578263963294167, -0.118573435441178,
    -0.165221384671569, -0.16900533418322, -0.168756839835035,
    -0.180391125223755, -0.180656844446289, -0.176773421796257,
    -0.185333909443916, -0.188228603437924, -0.1886657482609, -0.181768497164623,
    -0.171504941738718, -0.182633106078661, -0.178464686253981,
    -0.148940333601027, -0.197095206302847, -0.214655046781782,
    -0.205151106887526, -0.200993928570723, -0.181842959346968,
    -0.197221666529794, -0.141927125442055, -0.20814188074512,
    -0.176168014077669, -0.209113519059061, -0.173010232341285,
    -0.168374276986428, -0.181950989530223, -0.182572458277411,
    -0.19520507238331, -0.192910914684703, -0.187644986871977,
    -0.184732339747445, -0.193097075672569, -0.21368063379687,
    -0.149558615719411, -0.188506091050778, -0.149180605445445,
    -0.145144977091439, -0.139250327857729, -0.0819101284104886,
    -0.0898071362628176, 0.0737022220444434, 0.993688421218617,
    0.389542681708138, 0.259953546289957, 0.205350599725703, 0.0588315921898344,
    0.068310398426278, 0.222362991198682, 4.24720589521723, 3.21616278362766,
    3.98808332415544, 2.97829718431155, 3.57503688188657, 2.25378355544251,
    2.34552361719363, 2.48556086959792, 2.14164991010122, 2.41029827276356,
    1.72139977493018, 2.0652681668492, 2.53918511925901, 2.9049155500075,
    1.93596154510176, 2.69297446178442, 0.560338466617465, 1.40739127325152,
    1.84977483736185, 0.0311639670693444, 0.9244008754986, 2.54073692014061,
    3.49589753625571, 2.1020520602971, 3.39463560348713, 3.81614847006756,
    2.67336534517933, 3.46512866817994, 2.86223273103076, 3.34696925660175,
    2.31830876762886, 4.8690948596589, 3.77608263193865, 2.81594167937841,
    2.85669335000305, 3.36409124181523, 3.58666167478315, 2.5551643816567,
    2.871560481003, 2.59499431527887, 3.0460660339132, 2.62945085948493,
    3.38453888039937, 2.17545488086463, 3.08531722290761, 3.47569095635649,
    3.45332450091286, 4.08179484487432, 5.18647329191939, 3.8275016757487,
    3.43544655234586, 4.30531649667246, 5.7577025384763, 2.26254483943039,
    3.57883055703089, 3.81946175413703, 3.75710664313559, 3.65597842585424,
    2.95562452870516, 2.96818723577851, 4.96063879999325, 3.0152295886781,
    3.25530054681281, 3.04318400349156, 4.00477820124871, 2.53257581139844,
    3.83654843639179, 3.97156209414913, 3.19843913033136, 3.12841277281432,
    3.74901699325855, 2.58452610481333, 4.07948367440201, 3.49477793060226,
    2.55265605814176, 3.65172531214794, 4.21113524426459, 3.88614870405905,
    2.42952401985065, 3.67412557923281, 3.15622150965099, 4.56234678776478,
    3.00693421891971, 4.82948775687818, 2.82093569482146, 4.50702887582584,
    2.59072111843533, 4.15315923595639, 3.07625097183622, 4.86475639081858,
    2.92591510051391, 4.30928689286671, 2.56952272118732, 4.38054687944193,
    2.69361945216824, 4.35763507031079, 3.69351177973501, 3.12803604080913,
    5.44170695100145, 29.4259244672724, 32.9828115942182, 35.4799265405775,
    45.2030673847771, 35.8718187300156, 37.0015004833416, 31.8593602117472,
    39.5816772294549, 42.4565000997961, 39.6569084046984, 40.7235056388616,
    37.1900486216219, 30.4474700933519, 33.356837113987, 37.8107389646341,
    37.3308087730883, 41.311490008978, 30.6079338289173, 37.3141001920503,
    39.6074894304185, 31.5458462356564, 43.6345867936082, 36.3958336775302,
    33.5441094374619, 40.4904129264493, 39.3936295381436, 46.3730502969807,
    38.4278550998956, 39.6639738190179, 36.3430021227096, 40.7019498421003,
    44.1464773958682, 38.4077591711244, 39.8246594895966, 44.2231072084932,
    40.9214074311105, 39.8642828054064, 47.0045716594454, 37.886607189757,
    37.4707292040383, 44.9022590913998, 45.6087751860247, 37.827907312182,
    41.2744121034685, 47.7671605041138, 50.6406757144986, 44.5457651754922,
    46.5938747964513, 44.7833035842408, 42.8176393550137, 46.0649609560007,
    52.0308403160298, 49.2395360978872, 48.808163954923, 63.9055572907374,
    55.2600410928002, 81.6049287714722, 0.000451963765186214,
    0.000298242158776355, 0.0480072873638614, 0.0843186236798628,
    0.0681371436239281, 0.069044851875333, 0.0211798759620347,
    0.0281897790630876, 0.0104055751033274, 0.0120824966670239,
    0.011707655919303, 0.00429538925482804, 0.026171720990157, 0.012308844385037,
    0.0269759722595066, 0.0232562387976906, 0.0313134692773148,
    0.079915240555244, 0.000302009860188988, 0.000373995893675216,
    0.000327121024324469, 0.000289511073375572, 0.00254122844912876,
    0.00050351593941605, 0.0381012518540362, 0.061207199567179,
    0.0164527747364926, 0.0311321943242712, 0.0213906685773222,
    0.0109196467462609, 0.0531344836423717, 0.0558265639060381,
    0.00682380657158572, 0.00465712488456276, 0.00397286779636406,
    0.00277662050395859, 0.0185202939305158, 0.0930561707819025,
    0.0790291357609392, 0.016731201181585, 0.00531782522778706,
    0.00494229398710584, 0.00337672619689862, 0.00308754227749219,
    0.00366276021157812, 0.00494463653760302, 0.00287598887609978,
    0.00383351166717793, 0.0047708925478916, 0.0038543196942024,
    0.0104721924423347, 0.00291451323524859, 0.00362819379631916,
    0.0179679455696358, 0.00905614598651628, 0.00247207760679584,
    0.00447361168857941, 0.00366908522268333, 0.0326197409819614,
    0.00770454850178376, 0.0272379241761882, 0.00681812903812488,
    0.01359037878911, 0.015943675265893, 0.00747149837566544, 0.0047706447669614,
    0.00507289717154236, 0.00409950874599091, 0.00496426398708139,
    0.0157884348369869, 0.00326685795495063, 0.00961924408738495,
    0.00234769293601975, 0.0235644675360894, 0.0157442716894968,
    0.0226484911504306, 0.0243018893605037, 0.00675305289149495,
    0.00557860118108638, 0.0355142327777438, 0.0379576899918023,
    0.133688244189985, 0.0748662864611712, 0.0806387730805141,
    0.0126163133873244, 0.00944625252915908, 0.01067452573584,
    0.00566763983147886, 0.0784156319015868, 0.37368742725843, 0.301569531571255,
    0.345737490730525, 0.301882192686616, 0.319369338621276, 0.274727037936045,
    0.245797464828187, 0.259508702648762, 0.223004336733064, 0.217301404292581,
    0.199736152531824, 0.229196105255999, 0.24618493714455, 0.249911721281037,
    0.222692645788243, 0.272733429834378, 0.284751226053592, 0.23383405111903,
    0.23463803208804, 0.286957165325828, 0.279660982735174, 0.245679078463525,
    0.259884708719513, 0.253168514314394, 0.291808996034131, 0.31405070878821,
    0.304811067744619, 0.281843778965829, 0.30384017824912, 0.327130190338145,
    0.280930825881394, 0.339199181220197, 0.399706470227796, 0.30716459576377,
    0.356677012890964, 0.29307805423248, 0.275892591666876, 0.276905331791285,
    0.309983536420755, 0.254115617292309, 0.264757676009701, 0.2301075548516,
    0.245482961945492, 0.211116848535059, 0.263165308299423, 0.306396622283708,
    0.303431079014623, 0.35039086094755, 0.3585144837367, 0.326612916144276,
    0.271517331334452, 0.323424891345396, 0.352394906363902, 0.305721641390469,
    0.321660255799539, 0.389517670787107, 0.322656349155048, 0.281069028506717,
    0.325542685211865, 0.288295296162407, 0.341056275811543, 0.315975959969865,
    0.311253037058465, 0.280501757818375, 0.292343681542383, 0.264497862157222,
    0.337357830607314, 0.35108095997525, 0.330060637328088, 0.275611313969829,
    0.313649988337328, 0.349845744591972, 0.351822127715359, 0.30839037173216,
    0.293881021925633, 0.285003105182242, 0.323532112545353, 0.348505264029185,
    0.35471828546455, 0.282543609283627, 0.319776198783703, 0.361239520771941,
    0.305515298881325, 0.354453461507944, 0.30465498886816, 0.289440370900381,
    0.324868961491135, 0.327572119222309, 0.40198419191435, 0.32789933756631,
    0.356871929379416, 0.328049242596585, 0.358191054520311, 0.366871162547708,
    0.366584611470105, 0.335108039693735, 0.323552743225618, 0.326038687015209,
    0.347466139180839, 0.961126736588472, 1.14366231040271, 1.22002176883422,
    1.33878132207478, 1.1895940903257, 1.13949072782448, 1.31016537859454,
    1.51401113356501, 1.35339227350947, 1.42480567174781, 1.33572726933006,
    1.24760266491817, 1.14344391486348, 1.11460129165861, 1.16546980232285,
    1.25132426111194, 1.32089282961562, 1.0716783306619, 0.9886176473597,
    1.29744664519679, 1.2788325475053, 1.48926915692509, 1.15931481408171,
    1.13603404061275, 1.38404636736082, 1.17519672702231, 1.36203438277619,
    1.0041325777909, 1.23663978558971, 1.12906027845684, 1.26907367446955,
    1.39576829974855, 1.29965638129136, 1.17107528733454, 1.19373897164272,
    1.37216646260883, 1.00574688666852, 1.39892816108325, 1.18926906559097,
    1.12884087463772, 1.17472153728063, 1.2631449062306, 1.12624000033988,
    1.33200797667723, 1.36078599516355, 1.40091331389154, 1.09696129864518,
    1.58418975959604, 1.06910812831771, 1.32084785283526, 1.31399941478668,
    1.32160785312962, 1.47330656166298, 1.5869741852351, 1.53823108411503,
    1.29277322446939, 1.48029190033515, -0.00629183099714414,
    -0.00637015989860001, 0.228981690419934, 0.297853271490099,
    0.157162657613959, 0.167592997530028, 0.0688433902401993, 0.134691423734401,
    0.0399367441469645, 0.032001433567292, 0.0469146405435159,
    0.00541957241986157, 0.0627778935020487, 0.0214195242160578,
    0.0766545954005051, 0.0386898889449785, 0.0799419429183619,
    0.333144015677462, -0.00654404246164109, -0.00599147536580158,
    -0.00637783914887491, -0.00655351381888492, 0.00917655583560273,
    -0.00670036051977763, 0.167457836365323, 0.18784543527699,
    0.0414524917731458, 0.081491821947725, 0.059837945047583, 0.0168660878143587,
    0.219203408353813, 0.185541077479251, 0.0147313062752861,
    0.00664650291963209, 0.00548203958209514, 0.00243695577247904,
    0.0822235381820409, 0.216321947967178, 0.23904601586249, 0.0419241602561924,
    0.0106870696475267, 0.00948210281420048, 0.00333118946848288,
    0.00378418410422565, 0.00395140537312089, 0.00805617171583406,
    0.00253907349694194, 0.00474007132485133, 0.00876719212849908,
    0.00852487616560649, 0.0272196368043562, 0.00250540154580992,
    0.00342401423024419, 0.0500342643831151, 0.0152975755979989,
    0.000147078680621471, 0.0115292424768167, 0.00494021952262069,
    0.102365762147528, 0.04025797333323, 0.0659712047681977, 0.0329286653658769,
    0.0430210164569287, 0.0435809401719294, 0.0142918606105749,
    0.0139932898146529, 0.0101503871318722, 0.00506848616642697,
    0.00742985588980383, 0.0530834535721327, 0.00683247161070955,
    0.021138613889784, 0.00329806783273701, 0.0781625155463153,
    0.0581706111275482, 0.0656786698207845, 0.117695487698074,
    0.0174852309293054, 0.0153372506944967, 0.101001368129948, 0.111645643791326,
    0.32506413741177, 0.282756680677786, 0.200138561761029, 0.0592419095324135,
    0.0246826949940333, 0.0393082127174862, 0.029305500708011, 0.335146191678229,
    1.20739757526277, 1.07389234139091, 1.15744985562893, 0.790310219750312,
    0.973247148872008, 0.685528422443349, 0.659459051026891, 0.731101775065087,
    0.689211378178034, 0.507546334641913, 0.673829715109635, 0.715912557313984,
    0.722488796395769, 0.676392563782687, 0.645085484231646, 0.622432058816816,
    0.646034444991563, 0.730559893414332, 0.665784595508922, 0.724618551527804,
    0.579105851556671, 0.787847233364567, 0.642659415276265, 0.713026055981458,
    0.888487468705997, 0.880140807232507, 1.15459861121345, 0.698769906810594,
    1.25892761177425, 0.771418120711317, 1.03999344578324, 1.20211237375283,
    1.42028910326855, 0.827276934797293, 1.16309354903211, 0.803568339897337,
    0.673522950667835, 0.718368758913598, 1.20073948795713, 0.7974396658997,
    0.809759742396596, 0.784950639754503, 0.606712276797047, 0.701284469587344,
    0.809686979870337, 0.901842406575436, 1.03187409360212, 1.41153188262201,
    1.15568501185089, 1.07156944052045, 0.778793470725181, 0.943946295763116,
    1.05226016518028, 0.881474261877539, 1.17604466443923, 1.18535616457383,
    0.911387540025066, 0.86320725973363, 1.04472347249, 0.81930260642944,
    1.28973223350064, 0.866290082689299, 1.22360601712387, 0.918118266487066,
    0.922827583823046, 0.596997360728677, 0.818097392443362, 0.927030381381168,
    0.864055795948946, 0.786262332765934, 0.834278499379297, 1.22895207039192,
    1.32574800461002, 1.04279241692398, 0.764717532846417, 0.875848260254945,
    0.819725718562606, 0.955957478631073, 1.0864633133071, 0.775622330825551,
    0.841870608191207, 1.0875591682286, 0.784780899229886, 1.12562624209475,
    0.928550529020416, 0.914118705201696, 1.00858473456777, 1.03570370975613,
    1.37110464633801, 1.18545234123721, 1.20066532743248, 0.739317844343791,
    0.998398419428887, 1.17052277466514, 1.75170529526875, 1.23548863495702,
    1.0072702069923, 1.28844064738407, 1.25620953160106, 2.91625126018473,
    3.53549205585241, 4.09791842332011, 5.76424898746972, 4.80335304099393,
    3.95165643675341, 4.97792517891908, 5.36789038148374, 5.20817447660104,
    6.22637163094773, 5.50991458312464, 3.86920046897104, 4.01481756340974,
    3.33499054850013, 3.81457611796773, 4.66282809278498, 5.29790217364968,
    3.66770179696712, 3.13049133457211, 4.62105754013429, 4.99348971517442,
    5.98426768012582, 4.57217468964117, 4.13415066142565, 5.83388455402483,
    4.75009219465341, 6.59143115079033, 3.38860696621447, 4.71020218794852,
    4.83471670292554, 4.79450708557076, 5.5859730818723, 5.05584770848452,
    4.43910607688624, 4.62886848772941, 5.56140551771438, 3.5466257374051,
    6.37750185944153, 5.72249912493654, 4.39251425242875, 4.63271907738505,
    5.26407904750555, 4.16608615696053, 5.10197164999328, 4.82966904105771,
    5.50698262106079, 3.36768615150087, 6.18500022745634, 3.78712157574395,
    5.56343820107495, 5.68842586463965, 4.72452565715151, 6.14008591540721,
    6.18510201510592, 6.30107040251751, 5.10462439541111, 5.94581552047139,
    -0.00916139529056792, -0.00765436557733001, -0.107784964019195,
    -0.231191200361407, -0.197906376160917, -0.167373774714909,
    -0.054259629703606, -0.0992961406026165, -0.0323200779093616,
    -0.0413913648812143, -0.0298774420046674, -0.0135563187891902,
    -0.0919224208550233, -0.0459368407160752, -0.0700414653128224,
    -0.129044991413656, -0.135034302822314, -0.32866804540842,
    -0.00792452459236226, -0.00814175362838432, -0.00798264253950121,
    -0.00782897214688638, -0.0103571320052214, -0.0101166854467861,
    -0.0133429884049722, -0.166770840636848, -0.0311985146904562,
    -0.0718993728581322, -0.0619367662757523, -0.0225065221084477,
    -0.106459844480849, -0.142365651988203, -0.0181613994807466,
    -0.0165103977194748, -0.0129180011836775, -0.00873746568518652,
    -0.0366531961032349, -0.162453771998671, -0.183811632533907,
    -0.0724389328117577, -0.0168685181396984, -0.0174270217872661,
    -0.00980042262394876, -0.0114885669020467, -0.0119237166977976,
    -0.0207791996498354, -0.0129898140762442, -0.0112772511313847,
    -0.016235953421945, -0.0132998051980215, -0.0437706757886549,
    -0.00932397047040579, -0.0114741140863264, -0.0580988652961805,
    -0.0320847570854815, -0.00937155957884672, -0.017064468717315,
    -0.0132832678917968, -0.119966946432493, -0.0233697229845903,
    -0.0689683246160603, -0.0175059915961801, -0.0561292995323811,
    -0.0558530616375367, -0.0196003603555336, -0.0126040111434228,
    -0.0132390249306927, -0.020167600494613, -0.0177135290233887,
    -0.0568745589104799, -0.0119151889760317, -0.0462893023728372,
    -0.00800880220675226, -0.126512219209509, -0.0481520565795955,
    -0.065370729590836, -0.0871572322049784, -0.0189692019027818,
    -0.0121146768000954, -0.105028606514792, -0.108865230728757,
    -0.438207944980376, -0.1366775426091, -0.320429207567674,
    -0.0179377178095237, -0.0387699584532729, -0.0272106110542428,
    -0.0150639729562013, -0.261336989431156, -0.388157945381647,
    -0.412049778885185, -0.320831897271779, -0.411143898826891,
    -0.401294156077545, -0.438240181413441, -0.300326540766866,
    -0.334409939283914, -0.26208966468389, -0.26830174608452, -0.239422363733178,
    -0.293454394935285, -0.301024916439974, -0.321027229132113,
    -0.289952469033682, -0.36034506903387, -0.52093932564142, -0.411640050152749,
    -0.363751719062876, -0.486692875095925, -0.403130416132966,
    -0.405256173638255, -0.323219659772126, -0.353485525496241,
    -0.336938205785826, -0.392770804683944, -0.307817634708043,
    -0.338781391854424, -0.356015132620753, -0.402579165515095,
    -0.372764124034334, -0.28136950672778, -0.370385705137627,
    -0.422012290923314, -0.336251870462077, -0.379812240376685,
    -0.371441051894207, -0.351165995349105, -0.444246078568762,
    -0.363736509447893, -0.323929460220793, -0.337773965431155,
    -0.309462626272302, -0.388021919681506, -0.307533333675733,
    -0.326228895600688, -0.345821415668332, -0.440493083942472,
    -0.402382385259884, -0.328546528687808, -0.32816062263709,
    -0.384062770229033, -0.431181349058963, -0.420461008852381,
    -0.379734554281923, -0.419921834478709, -0.430592393084476,
    -0.321480025450584, -0.352646394242016, -0.379273025784121,
    -0.320185361461885, -0.375068931742287, -0.285688019534665,
    -0.298996294988579, -0.31695133994016, -0.383709875033578,
    -0.500343640609472, -0.512764013279896, -0.405826898706391,
    -0.378997749607584, -0.43383828434994, -0.376190583177842,
    -0.388931504803331, -0.357399143648326, -0.319827392610684,
    -0.373058470590037, -0.380938980129992, -0.423198945521155,
    -0.487914749791299, -0.370173313369901, -0.398552800972669,
    -0.310207331105373, -0.393818494894227, -0.398929676850394,
    -0.377458606995538, -0.423851956578052, -0.433851757506572,
    -0.421678673876596, -0.460742574670196, -0.358599821509089,
    -0.348284243621751, -0.402663264725994, -0.468842839083677,
    -0.495203658205066, -0.392291204734205, -0.396253858065874,
    -0.333310494520421, -0.36667790799295, -0.347511589845891,
    -0.770837234011864, -0.895877100935804, -0.879081524649439,
    -0.88607106836547, -0.700248082298628, -0.878863933732385,
    -0.707808816542336, -0.828281759285298, -0.792688688931428,
    -0.889391272853463, -0.711680638593714, -0.924466753268476,
    -0.833221035589461, -0.816913820353629, -0.603863175847747,
    -0.780349542358818, -0.832451008419022, -0.954903392272828,
    -0.888988432156231, -0.734074476728476, -0.784686642405841,
    -0.854868591135917, -0.52130692545968, -0.89328998673569, -0.76808469842316,
    -0.628891013696258, -0.702776863304914, -0.907393763199671,
    -0.803314455281585, -0.746236187961955, -0.893230346988993,
    -0.839047086646552, -0.926251400377845, -0.734745315067091,
    -0.751982519763878, -0.863203081692536, -0.692273786410333,
    -0.633615867741087, -0.463445855843162, -0.749805275252324,
    -0.636481275441413, -0.521528260002451, -0.599819701615391,
    -0.722316332150902, -0.812346967811988, -0.581341664561544,
    -0.868658914364143, -0.632416026376594, -0.625294924013276,
    -0.547940573040295, -0.318081041197452, -0.799311076121852,
    -0.87066569149278, -0.631977585262865, -0.767710724752195,
    -0.726491406092615, -0.650529781259422, -0.342405809251406,
    -0.337279675945309, -0.106327810697275, -0.158224598521657,
    -0.0271529983005288, 0.293106079148189, -0.134291969790532,
    -0.0756813410508732, -0.170786998167568, -0.196076932710698,
    -0.181356411573579, -0.192466355110217, -0.175943933064245,
    -0.184629512722152, -0.253846385107435, -0.203235092430561,
    -0.188434315537133, -0.532140092625756, -0.351603805855075,
    -0.345956761746581, -0.351190868090352, -0.351038761343783,
    -0.330661830194012, -0.36241671090587, 0.278678978749028, -0.539711498570492,
    -0.112934279182954, -0.117988682341662, -0.161095598103117,
    -0.161025819798973, -0.0991946485854701, -0.122626963263036,
    -0.158745242779326, -0.172235273038471, -0.171947937986705,
    -0.169027850241014, -0.140926479898964, 0.357838076412855, 0.608893965019851,
    -0.123773193623336, -0.16561048596812, -0.164382634329419,
    -0.166345071843829, -0.177868893898343, -0.173317474797937,
    -0.176022229372947, -0.182997317836965, -0.180798277881699,
    -0.187749019562661, -0.180890121874806, -0.175949236932859,
    -0.179736774075305, -0.174823681816529, -0.152198939773389,
    -0.18562596534214, -0.210944052231684, -0.199106388110678,
    -0.195579770605385, -0.18062466086518, -0.193180956747049,
    -0.119427976730769, -0.203509472238879, -0.164309971010146,
    -0.202998502716338, -0.167875631899814, -0.163112261294292,
    -0.176417562478208, -0.178349192397569, -0.191759662581594,
    -0.19363261192769, -0.182769228531672, -0.181721712049437,
    -0.190493723485731, -0.194665654159812, -0.162819418053741,
    -0.184732733119413, -0.145284797688599, -0.14271994334288,
    -0.144566986345344, -0.0839179260646422, -0.0984001304337133,
    -0.0704094542231027, 0.920648852184275, 0.490901810285805, 0.238993060029918,
    0.189340239281242, 0.0654098688409682, 0.0668757990853168, 0.211777359840398,
    4.11715279411116, 3.1886722609799, 4.01708422669355, 2.66709321519585,
    3.74799071584902, 2.11901395748789, 1.97897445136316, 2.33272600079106,
    1.90293780980691, 2.11444595229815, 1.73981069858459, 1.89074965502213,
    2.31097299858721, 2.97286103897225, 1.66870372239606, 2.56705359511122,
    0.844541484349054, 1.64166802886222, 1.8097026280157, 0.00781207032541342,
    1.18700994897948, 2.6565553098344, 3.66481803056936, 2.39343599709227,
    3.23276578644064, 3.8615445459342, 2.85731740980108, 3.23372824734535,
    2.90595596954844, 3.16055250326269, 2.22355227904068, 4.33815140059522,
    3.47254938810895, 3.08203376082199, 2.93935299382011, 3.05250045779209,
    3.51662399889908, 2.4942310005602, 3.07125704190881, 2.69864727986873,
    3.15493520008894, 2.80631757279575, 3.50696563499646, 2.03061465595928,
    3.00315378983383, 3.72750822662089, 3.59079574891748, 3.83886025649069,
    4.39466474621014, 3.76239600151466, 3.71587151391252, 4.44349404126688,
    5.72807191923134, 2.56139654611775, 3.7504012131588, 3.22241743929892,
    3.40869905420392, 3.2639909205019, 2.9283319911514, 2.71636295501483,
    4.30187925279305, 2.97197013626913, 3.02430995599471, 3.00841865092493,
    4.00853818863823, 2.80762537110052, 3.96789285092546, 3.83013806983957,
    2.95349421688033, 2.9186110031591, 3.5403177988428, 2.75744077650249,
    3.72278369227332, 3.40335558374541, 2.71812965037006, 3.53834167852093,
    3.92786722563706, 3.49558231605565, 2.4078855010167, 3.55711094820419,
    2.88938228608127, 4.61281566805482, 2.95629344032919, 4.61106650028523,
    2.87201360673364, 4.43377108973695, 2.49905636648777, 4.00737258777281,
    3.4044234861475, 4.5638879891862, 3.02542457611229, 3.76667579628371,
    2.75955827778939, 4.10581110353222, 2.98781627230194, 4.04559972134318,
    3.56747229226355, 2.90842421300586, 5.01753203108475, 28.656316750191,
    32.8273245460202, 34.5824308104601, 43.5860779301054, 34.4729373976474,
    35.6119246160946, 31.8519857444026, 39.6903341917286, 41.9418392101259,
    38.9397744959233, 39.2983176694736, 35.5592325979822, 30.4893344554993,
    33.4421791287291, 37.5837242531644, 36.1889279383811, 38.2818555649409,
    29.9656944073967, 36.4239099750158, 39.1641178816357, 30.9429713160882,
    42.7640495551054, 34.8743769637091, 32.0890233971647, 38.6526169787827,
    37.7849480249902, 44.7705216053421, 36.7688269288749, 38.4862247878429,
    35.6084345135905, 40.2097936362177, 42.32980387737, 37.9371929225819,
    38.302492471624, 43.1433740830749, 39.6802296999627, 39.9357377184739,
    44.1326286635952, 37.1050671531373, 37.1143857151029, 43.4847440339207,
    44.3008835452229, 36.6824924013357, 39.9829305165219, 45.8272105622544,
    50.3926508412704, 43.1455615393743, 45.4885531636656, 44.3439129261837,
    40.7493114796369, 45.8653667185415, 50.9969815124154, 49.3863365528443,
    48.120825391347, 62.9762928248327, 54.866138790752, 80.0075122360497, 0.0,
    0.0, 1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 13.0, 8.0,
    10.0, 9.0, 10.0, 9.0, 6.0, 12.0, 8.0, 10.0, 7.0, 7.0, 7.0, 7.0, 9.0, 9.0,
    9.0, 7.0, 7.0, 8.0, 8.0, 9.0, 9.0, 10.0, 11.0, 13.0, 9.0, 9.0, 10.0, 10.0,
    6.0, 10.0, 10.0, 9.0, 7.0, 10.0, 11.0, 9.0, 7.0, 9.0, 12.0, 7.0, 9.0, 5.0,
    7.0, 10.0, 10.0, 9.0, 11.0, 11.0, 10.0, 10.0, 13.0, 9.0, 9.0, 9.0, 11.0, 9.0,
    9.0, 8.0, 9.0, 10.0, 10.0, 13.0, 11.0, 8.0, 10.0, 10.0, 11.0, 10.0, 8.0,
    10.0, 9.0, 11.0, 9.0, 8.0, 10.0, 10.0, 13.0, 9.0, 7.0, 13.0, 9.0, 11.0, 11.0,
    11.0, 10.0, 13.0, 9.0, 10.0, 11.0, 12.0, 12.0, 11.0, 7.0, 10.0, 9.0, 12.0,
    10.0, 10.0, 9.0, 13.0, 14.0, 10.0, 15.0, 10.0, 11.0, 12.0, 13.0, 13.0, 14.0,
    11.0, 14.0, 10.0, 14.0, 12.0, 13.0, 14.0, 10.0, 11.0, 14.0, 12.0, 12.0, 11.0,
    13.0, 16.0, 13.0, 15.0, 11.0, 14.0, 14.0, 14.0, 14.0, 13.0, 14.0, 16.0, 16.0,
    16.0, 15.0, 18.0, 11.0, 13.0, 14.0, 14.0, 16.0, 16.0, 16.0, 14.0, 12.0, 16.0,
    15.0, 16.0, 17.0, 14.0, 16.0, 17.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 30.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 46.0, 37.0, 46.0, 37.0, 43.0, 36.0, 38.0, 41.0,
    38.0, 40.0, 42.0, 40.0, 41.0, 47.0, 33.0, 37.0, 41.0, 45.0, 36.0, 38.0, 42.0,
    44.0, 42.0, 42.0, 39.0, 46.0, 40.0, 43.0, 47.0, 46.0, 28.0, 40.0, 41.0, 41.0,
    38.0, 38.0, 46.0, 39.0, 39.0, 40.0, 42.0, 41.0, 42.0, 40.0, 39.0, 39.0, 43.0,
    40.0, 46.0, 43.0, 42.0, 42.0, 43.0, 40.0, 47.0, 38.0, 43.0, 45.0, 38.0, 38.0,
    45.0, 39.0, 38.0, 47.0, 47.0, 41.0, 42.0, 44.0, 40.0, 47.0, 39.0, 39.0, 46.0,
    39.0, 40.0, 45.0, 37.0, 43.0, 41.0, 36.0, 38.0, 46.0, 37.0, 46.0, 39.0, 47.0,
    37.0, 47.0, 37.0, 45.0, 38.0, 44.0, 40.0, 45.0, 40.0, 45.0, 40.0, 40.0, 39.0,
    41.0, 46.0, 42.0, 45.0, 41.0, 44.0, 44.0, 44.0, 45.0, 41.0, 44.0, 45.0, 41.0,
    47.0, 46.0, 43.0, 40.0, 38.0, 45.0, 46.0, 38.0, 46.0, 42.0, 37.0, 37.0, 44.0,
    46.0, 40.0, 45.0, 42.0, 46.0, 43.0, 47.0, 43.0, 45.0, 45.0, 44.0, 45.0, 45.0,
    41.0, 43.0, 41.0, 42.0, 40.0, 45.0, 45.0, 45.0, 41.0, 40.0, 41.0, 43.0, 46.0,
    44.0, 46.0, 45.0, 46.0, 46.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.5, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.5, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 3.53846153846154, 4.625, 4.6, 4.11111111111111, 4.3,
    4.0, 6.33333333333333, 3.41666666666667, 4.75, 4.0, 6.0, 5.71428571428571,
    5.85714285714286, 6.71428571428571, 3.66666666666667, 4.11111111111111,
    4.55555555555556, 6.42857142857143, 5.14285714285714, 4.75, 5.25,
    4.88888888888889, 4.66666666666667, 4.2, 3.54545454545455, 3.53846153846154,
    4.44444444444445, 4.77777777777778, 4.7, 4.6, 4.66666666666667, 4.0, 4.1,
    4.55555555555556, 5.42857142857143, 3.8, 4.18181818181818, 4.33333333333333,
    5.57142857142857, 4.44444444444445, 3.5, 5.85714285714286, 4.66666666666667,
    8.0, 5.57142857142857, 3.9, 4.3, 4.44444444444445, 4.18181818181818,
    3.90909090909091, 4.2, 4.2, 3.30769230769231, 4.44444444444445,
    5.22222222222222, 4.22222222222222, 3.90909090909091, 5.0, 4.22222222222222,
    4.75, 5.0, 3.9, 3.8, 3.61538461538462, 4.27272727272727, 5.125, 4.2, 4.4,
    3.63636363636364, 4.7, 4.875, 3.9, 5.11111111111111, 3.54545454545455,
    4.44444444444445, 5.625, 3.7, 4.3, 3.15384615384615, 4.0, 5.42857142857143,
    3.53846153846154, 4.11111111111111, 4.18181818181818, 3.54545454545455,
    4.27272727272727, 3.7, 3.61538461538462, 4.11111111111111, 4.5,
    3.45454545454545, 3.66666666666667, 3.33333333333333, 4.09090909090909,
    5.71428571428571, 4.5, 4.44444444444445, 3.33333333333333, 3.9, 4.1,
    5.11111111111111, 3.23076923076923, 3.21428571428571, 4.1, 2.93333333333333,
    4.4, 4.0, 3.75, 3.15384615384615, 3.38461538461538, 3.21428571428571,
    3.72727272727273, 3.35714285714286, 4.6, 3.07142857142857, 3.33333333333333,
    2.92307692307692, 3.21428571428571, 4.6, 3.45454545454545, 3.28571428571429,
    3.5, 3.08333333333333, 3.36363636363636, 3.38461538461538, 2.875,
    3.07692307692308, 3.0, 3.81818181818182, 3.28571428571429, 3.07142857142857,
    3.35714285714286, 3.07142857142857, 3.46153846153846, 3.21428571428571, 2.75,
    2.8125, 2.8125, 2.73333333333333, 2.38888888888889, 3.72727272727273,
    3.23076923076923, 2.85714285714286, 3.21428571428571, 2.8125, 2.8125, 2.5625,
    2.85714285714286, 3.41666666666667, 2.6875, 3.06666666666667, 2.75,
    2.70588235294118, 3.21428571428571, 2.875, 2.70588235294118, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 2.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    1.0, 1.0, 0.0, 1.0, 1.0, 14.0, 15.0, 17.0, 19.0, 14.0, 15.0, 14.0, 15.0,
    16.0, 15.0, 17.0, 16.0, 13.0, 14.0, 16.0, 16.0, 19.0, 13.0, 17.0, 17.0, 14.0,
    16.0, 15.0, 16.0, 18.0, 17.0, 17.0, 18.0, 18.0, 14.0, 16.0, 18.0, 14.0, 16.0,
    14.0, 15.0, 15.0, 18.0, 13.0, 13.0, 16.0, 16.0, 16.0, 16.0, 17.0, 18.0, 18.0,
    15.0, 19.0, 16.0, 15.0, 18.0, 15.0, 16.0, 21.0, 19.0, 31.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 7.0, 8.0, 7.0,
    12.0, 4.0, 5.0, 5.0, 4.0, 4.0, 2.0, 4.0, 7.0, 5.0, 4.0, 6.0, 4.0, 3.0, 3.0,
    4.0, 6.0, 5.0, 8.0, 6.0, 7.0, 7.0, 6.0, 8.0, 5.0, 7.0, 5.0, 6.0, 3.0, 6.0,
    6.0, 6.0, 9.0, 6.0, 4.0, 4.0, 6.0, 4.0, 6.0, 5.0, 8.0, 8.0, 7.0, 5.0, 11.0,
    6.0, 6.0, 9.0, 9.0, 6.0, 7.0, 9.0, 8.0, 8.0, 9.0, 8.0, 6.0, 7.0, 6.0, 5.0,
    7.0, 5.0, 8.0, 9.0, 7.0, 8.0, 7.0, 4.0, 7.0, 8.0, 7.0, 9.0, 9.0, 7.0, 6.0,
    6.0, 8.0, 8.0, 6.0, 9.0, 4.0, 8.0, 6.0, 7.0, 6.0, 7.0, 6.0, 9.0, 7.0, 9.0,
    6.0, 7.0, 6.0, 7.0, 6.0, 9.0, 7.0, 4.0, 7.0, 9.0, 12.0, 7.0, 8.0, 12.0, 9.0,
    7.0, 11.0, 8.0, 11.0, 5.0, 8.0, 6.0, 10.0, 10.0, 4.0, 9.0, 9.0, 9.0, 5.0,
    6.0, 7.0, 9.0, 9.0, 6.0, 14.0, 10.0, 9.0, 10.0, 11.0, 13.0, 9.0, 15.0, 12.0,
    13.0, 15.0, 15.0, 6.0, 7.0, 9.0, 11.0, 9.0, 13.0, 11.0, 11.0, 8.0, 10.0,
    13.0, 13.0, 11.0, 14.0, 16.0, 10.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 41.0, 43.0, 42.0, 43.0, 38.0, 46.0, 45.0, 45.0, 46.0,
    46.0, 48.0, 46.0, 43.0, 45.0, 46.0, 44.0, 46.0, 47.0, 47.0, 46.0, 44.0, 45.0,
    42.0, 44.0, 43.0, 43.0, 44.0, 42.0, 44.0, 43.0, 45.0, 43.0, 45.0, 44.0, 44.0,
    44.0, 41.0, 44.0, 45.0, 46.0, 44.0, 46.0, 44.0, 45.0, 42.0, 42.0, 43.0, 44.0,
    39.0, 44.0, 44.0, 41.0, 41.0, 44.0, 43.0, 41.0, 42.0, 42.0, 41.0, 42.0, 43.0,
    43.0, 43.0, 45.0, 43.0, 45.0, 42.0, 41.0, 43.0, 42.0, 43.0, 45.0, 42.0, 42.0,
    43.0, 41.0, 41.0, 43.0, 44.0, 44.0, 42.0, 42.0, 44.0, 41.0, 46.0, 42.0, 44.0,
    43.0, 42.0, 43.0, 43.0, 41.0, 43.0, 41.0, 43.0, 42.0, 44.0, 42.0, 43.0, 27.0,
    28.0, 29.0, 24.0, 27.0, 23.0, 29.0, 27.0, 22.0, 26.0, 26.0, 23.0, 29.0, 25.0,
    29.0, 26.0, 25.0, 27.0, 23.0, 29.0, 27.0, 25.0, 26.0, 29.0, 26.0, 26.0, 24.0,
    23.0, 26.0, 22.0, 24.0, 23.0, 26.0, 23.0, 23.0, 26.0, 20.0, 20.0, 24.0, 22.0,
    19.0, 28.0, 27.0, 25.0, 22.0, 23.0, 19.0, 24.0, 20.0, 26.0, 25.0, 19.0, 22.0,
    23.0, 15.0, 15.0, 9.0, 0.0, 0.000874264526367186, 13.8833206787109,
    75.2715529266358, 5.90390834655761, 16.940623727417, 5.15641217651367,
    1.30702546691894, 2.12271427001953, 0.629470458984377, 3.45334487915039,
    1.63400039978027, 5.41169741821289, 2.6831178314209, 2.72945385131836,
    0.210697750854491, 4.05921019592285, 16.5891693878174, 0.0227308776855469,
    0.0, 0.0236051422119141, 0.0236051422119141, 0.0935463043212891,
    0.0236051422119141, 38.9895750823975, 21.7936661132812, 2.12271427001953,
    2.56771491394043, 4.83031150817871, 1.16539461364746, 20.8136155792236,
    5.46153049621582, 2.12271427001953, 0.675806478881836, 1.18987402038574,
    0.117151446533203, 6.27547077026367, 11.3173542938233, 105.069985043335,
    6.44070676574707, 1.93649592590332, 0.348831546020508, 0.280638912963867,
    0.326974932861328, 0.163487466430664, 0.536798419189453, 0.210697750854492,
    0.0699411621093751, 0.280638912963867, 0.373310952758789, 1.70219303283691,
    0.0699411621093751, 1.21347916259766, 1.86655476379394, 3.45159635009766,
    0.093546304321289, 1.35336148681641, 3.77944554748535, 14.5827322998047,
    7.34994187316895, 2.09998339233399, 2.61405093383789, 1.02726081848145,
    0.559529296874999, 1.79661360168457, 0.0708154266357424, 0.373310952758789,
    0.886504229736328, 1.1671431427002, 0.583134439086914, 0.956445391845703,
    0.420521237182617, 0.862899087524414, 20.8354721923828, 4.13002562255859,
    10.5707323883057, 18.6200858825684, 1.07359683837891, 1.75027758178711,
    0.489588134765626, 0.232554364013673, 20.9771030456543, 19.2032203216553,
    8.39993356933594, 3.89572272949219, 2.93927733764648, 4.45612629089355,
    4.03560505371094, 22.983540133667, 0.0463360198974669, 6.46343764343264,
    1.67946215515136, 27.0663554718018, 1.70306729736328, 1.86917755737304,
    7.53703448181153, 0.0463360198974669, 2.63678181152345, 4.1300256225586,
    0.513193276977539, 8.702429095459, 6.88308461608887, 1.68208494873048,
    10.3145728820801, 13.2538502197266, 19.8799010650635, 4.96931956787109,
    4.76037034606934, 29.9837761962891, 34.6506002380371, 8.77324452209474,
    25.5259013763428, 15.2130770233154, 5.31989964294434, 5.90303408203124,
    7.8412785369873, 12.7869929626465, 23.0762121734619, 31.0328936279297,
    31.8966669799805, 24.7338177154541, 42.3030376373291, 25.7365991271973,
    11.8760093261719, 4.9693195678711, 20.9290184967041, 13.8142537811279,
    1.77388272399902, 25.3388087677002, 3.26625227050781, 13.3456479949951,
    5.71594147338868, 11.6661858398438, 9.14568121032714, 9.56707671203613,
    6.36989133911133, 7.11476471557616, 8.32999240722656, 4.71403432617187,
    8.51795928039551, 39.7353227233887, 6.93116916503908, 12.4818746429443,
    8.3772026916504, 10.7805558746338, 0.489588134765626, 17.1268420715332,
    12.5990260894775, 18.5737498626709, 25.2234058502197, 2.40335318298339,
    6.16006785278321, 12.6243797607422, 4.41066453552246, 7.21180807800293,
    5.29629450073243, 11.4790932312012, 0.604991052246092, 6.97575665588379,
    15.469236529541, 17.3602707000732, 23.2405739044189, 12.4836231719971,
    19.8090856384277, 8.72778276672364, 30.9402215881348, 4.8303115081787,
    5.27268935852052, 15.3529593475342, 12.3892026031494, 16.4964973480225,
    11.9232196105957, 12.5063540496826, 29.2808675170898, 4.71316006164551,
    33.9258349456787, 5.24995848083498, 10.3136986175537, 4.78310122375488,
    1.19074828491212, 8.00389173889161, 17.8498588348389, 1.77213419494628,
    4.97019383239746, 10.7106147125244, 10.0575391113281, 13.8605898010254,
    3.66404263000489, 54.0855006591797, 47.4139880584717, 84.1619488952637,
    124.575700891113, 31.1500450744629, 141.726148104858, 102.479539251709,
    90.0186469573974, 71.1205449554443, 56.4197869445801, 39.5727095214844,
    15.5409262207031, 78.3297302398682, 40.2966005493164, 117.786162579346,
    98.1860261627197, 61.5071322235107, 58.3090725860596, 64.2593169525146,
    72.0297800628663, 108.220834396362, 74.1061583129883, 93.6835638519287,
    90.0903366485596, 62.3000901489258, 87.0575130065918, 75.1797551513672,
    67.712661831665, 117.342036199951, 117.366515606689, 73.4766878540039,
    119.139524066162, 12.3664717254639, 41.602751751709, 86.588907220459,
    132.252617697144, 16.3338841461182, 106.819388360596, 33.9039783325195,
    22.7728423828125, 9.80050534057619, 40.552760055542, 1.89015990600586,
    42.162281048584, 47.8790967864991, 27.6958259307861, 11.457236618042,
    33.4135159332275, 4.29351308898924, 99.1896818389892, 25.4804396209717,
    22.7256320983887, 81.572377368164, 69.3693931091309, 64.7733844940185,
    54.086374923706, 28.5814558959961, 0.0, -3.80115011463994E-5,
    13.8833206787109, 25.0905176422119, 0.236156333862305, -0.705859321975708,
    5.15641217651367, -0.435675155639648, 1.06135713500976, -0.0251788183593751,
    1.7266724395752, -0.0817000199890137, -0.901949569702149, -2.6831178314209,
    -2.72945385131836, 0.00540250643216644, 1.01480254898071, 16.5891693878174,
    0.00174852905273437, 0.0, 0.0236051422119141, 0.000907890085073617,
    -0.0935463043212891, -0.0236051422119141, 6.49826251373291,
    -10.8968330566406, 0.707571423339844, 0.213976242828369, 1.61010383605957,
    -0.58269730682373, -0.867233982467651, 2.73076524810791, -2.12271427001953,
    0.0375448043823242, -1.18987402038574, -0.0146439308166504,
    -0.348637265014648, 11.3173542938233, 52.5349925216675, -3.22035338287354,
    -1.93649592590332, -0.0174415773010254, 0.093546304321289,
    -0.326974932861328, 0.163487466430664, -0.536798419189453, 0.006384780328924,
    -0.0699411621093751, 0.280638912963867, -0.373310952758789,
    -0.851096516418457, 0.0699411621093751, -0.0275790718772195,
    1.86655476379394, 1.15053211669922, 0.093546304321289, -1.35336148681641,
    -1.25981518249512, -2.91654645996094, -0.204165052032471, 2.09998339233399,
    1.30702546691895, -0.256815204620361, -0.0932548828124998, 1.79661360168457,
    -0.0236051422119141, -0.12443698425293, -0.886504229736328, -1.1671431427002,
    -0.583134439086914, 0.956445391845703, 0.420521237182617, 0.215724771881104,
    -20.8354721923828, -4.13002562255859, -0.310903893773696, -9.31004294128418,
    -1.07359683837891, -1.75027758178711, 0.023313720703125, -0.0387590606689455,
    10.4885515228272, -0.800134180068969, 0.399996836635045, -0.0905982030114462,
    2.93927733764648, 2.22806314544678, 2.01780252685547, 22.983540133667,
    0.0463360198974669, 0.239386579386394, -0.559820718383788, -9.02211849060059,
    -0.0630765665690105, 0.143782889028695, -0.203703634643555,
    0.00272564822926276, -0.439463635253908, 0.137667520751953,
    0.256596638488769, -1.45040484924317, -0.688308461608887, -0.168208494873048,
    0.606739581298829, 0.441795007324219, -0.946661955479214, -0.828219927978515,
    0.297523146629334, 4.99729603271484, 5.77510003967285, -0.381445414004119,
    -8.50863379211426, 1.52130770233154, 5.31989964294434, 0.13415986550071,
    0.280045662035261, -1.42077699584961, -7.6920707244873, 0.862024822998047,
    -1.2267948838454, 0.727465226925121, 4.70033751525879, 1.83832850908552,
    1.31955659179687, 0.552146618652344, 0.597971957048689, 1.38142537811279,
    -0.295647120666504, -0.589274622504656, -0.296932024591619, 1.33456479949951,
    -0.571594147338868, -1.94436430664063, 0.914568121032714, 1.59451278533936,
    0.193027010276101, -0.296448529815673, 0.640768646709735, -1.57134477539062,
    -1.41965988006592, 2.83823733738491, -0.330055674525671, -4.16062488098145,
    -2.7924008972168, 1.34756948432922, 0.037660625751202, 1.22334586225237,
    -2.09983768157959, -0.640474133195548, -1.94026198847844, -0.801117727661132,
    -0.280003084217419, 0.971106135441707, 4.41066453552246, 0.515129148428781,
    0.331018406295777, -0.546623487200056, 0.0672212280273436,
    -0.317079847994718, -1.5469236529541, 0.642972988891602, 2.58228598937988,
    12.4836231719971, -1.98090856384277, 0.545486422920227, -2.38001704524114,
    -0.178900426228841, -1.0545378717041, 1.91911991844177, -1.3765780670166,
    16.4964973480225, 1.32480217895508, 0.833756936645507, -14.6404337585449,
    -0.10246000134012, -3.76953721618652, 0.583328720092775, 1.14596651306152,
    0.18396543168288, -0.119074828491212, 0.88932130432129, 2.23123235435486,
    0.136318014995868, -0.177506922585623, -0.446275613021851, 0.773656854717548,
    1.06619921546349, 0.0796531006522801, -54.0855006591797, 5.26822089538574,
    3.8255431316029, -124.575700891113, -7.78751126861572, -4.88710855533994,
    4.65816087507768, -90.0186469573974, -2.45243258467049, -18.8065956481934,
    1.8844147391183, 1.19545586313101, 7.83297302398682, -2.51853753433228,
    5.35391648087935, -98.1860261627197, 6.8341258026123, 6.47878584289551,
    -64.2593169525146, 8.00330889587403, 12.0245371551514, 3.36846174149947,
    10.4092848724365, -5.63064604053497, -20.7666967163086, 2.48735751447405,
    -4.69873469696045, 3.22441246817453, -117.342036199951, 13.0407239562988,
    -4.59229299087524, -29.7848810165405, -12.3664717254639, 1.98108341674805,
    2.62390627940785, -4.72330777489798, -5.44462804870605, 5.34096941802978,
    4.23799729156494, 2.5303158203125, 0.296985010320491, 4.50586222839356,
    -0.630053302001954, 4.2162281048584, 1.450881720803, 9.2319419769287,
    1.27302629089356, 4.17668949165344, -0.159019003295898, 3.09967755746841,
    3.18505495262146, 1.13628160491943, -3.02119916178385, 7.70771034545899,
    -1.85066812840053, 1.80287916412353, 1.02076628199986, 0.429263882446289,
    0.429263882446289, 14.0581735839844, 578.499088568115, 36.2164080047607,
    14.3108360321045, -94.7370526062012, -7.79581678161621, -64.3834625152588,
    5.84533262329102, 2.51525904235839, -22.7142666595459, 29.7249938964844,
    -34.5177120300293, 16.0016636260986, -8.1743733215332, 34.7187928710937,
    -376.31405140686, 1.0578600769043, 0.257908035278321, -0.211572015380859,
    -0.490462399291992, -0.143379382324219, 0.139008059692383, 346.406336224365,
    377.031822583008, -16.1161922790527, -101.719803378296, 40.1156277923584,
    -30.9795634918213, 1.88404005432129, 48.56714296875, -30.2338158508301,
    -7.09640516052246, 2.42695832519531, -3.07915966186523, -71.0191302703857,
    925.922194436646, -583.425569174194, -3.97353227233887, -24.4304479248047,
    -3.4052603302002, 3.82490730285645, 1.67683936157227, -0.234302893066406,
    0.0227308776855466, -13.6490177856445, -1.21610195617676, -2.31155540771484,
    0.766729989624024, 1.8945312286377, -20.5137428466797, -6.18717005310058,
    -70.0163488586426, -10.8627367401123, 1.70394156188965, -1.50286072082519,
    3.84588965148926, 35.4418096343994, 8.18311596679688, -0.643458691406244,
    -26.090676260376, 20.6920928100586, -84.3376760650635, -11.9922865081787,
    -15.9439621673584, 9.81799063110351, -8.84930553588867, -5.48950696105957,
    21.2542449005127, -12.9321208740234, -6.18804431762696, -3.34843313598633,
    -4.25067412719726, 0.537672683715832, 40.8115423553467, 131.529600933838,
    -2.2433627746582, 5.33825919799804, 77.7929318206787, 18.2730028656006,
    -235.799633935547, -49.2403266540527, 91.4786687164307, -9.56969950561522,
    -2.98998468017578, 0.975679211425779, 10.1913015838623, -77.757961239624,
    252.556662112427, 404.673444113159, -22.1153954589844, -13.4829075256348,
    160.534200860596, 229.574870507812, -216.849950326538, 217.050156903076,
    -130.285522512817, 300.018734719849, 208.509466744995, 414.097141442871,
    -157.980474179077, 267.402548034668, 101.181256430054, 345.287277630615,
    235.928150820923, 221.458198645019, 381.476583435059, 136.590718276978,
    231.027023886108, -239.207517059326, 83.0026741333007, 793.909125219727,
    279.367732342529, 15.2629101013184, -104.132773471069, 255.384907855225,
    -275.113561157227, 416.12368661499, -191.324923214722, 323.439407116699,
    392.814920077515, -55.7378606140137, 181.487698764038, -190.530216760254,
    200.266900790405, -16.4720179412842, -36.4419682525635, 200.645457330322,
    -86.1386609893798, 411.294249371338, 38.012147341919, 523.529706472778,
    597.913880905151, 78.7502514770508, 414.088398797607, 203.634567745972,
    -165.845357858276, 102.318674578857, 186.998188073731, -119.344976229858,
    190.30727930603, 171.276289096069, -189.87888968811, 31.3152810699463,
    316.77051730957, -335.354758346558, 15.0076248596192, 182.698555133057,
    -89.3603257690428, -35.7679103027344, 254.610309484863, 26.0670711181641,
    391.572590185547, 162.919194488525, -202.366009918213, 153.65024197998,
    246.261083258057, -135.746178744507, 538.292537265015, 174.209446582031,
    -119.022372619629, 107.73037199707, 182.749262475586, -210.623438369751,
    643.220017190552, 201.297658666992, 77.1381076904297, -116.708194418335,
    60.2228376342774, 163.26715177002, 162.168201260376, -43.9588946502686,
    -51.0483056945801, 178.053587704468, 175.095950811768, -62.2275261932373,
    -96.1979486297607, 146.721695608521, 160.931116955566, -48.5225554779053,
    -1.86742902832027, 250.015175134278, 183.025530065918, -28.9967315460204,
    98.5575885864257, 247.356536709595, -245.548557669067, 70.550524484253,
    720.823233609009, 157.608911755371, 70.9002302947998, 70.5050627288819,
    32.3818837921143, 191.957016467285, 187.122333636475, 132.923178588867,
    -155.313093109131, 87.4727886566163, 373.791798248291, 390.17114414978,
    134.836069372559, 62.9225664916991, 264.498241278076, 34.5727906951904,
    133.436371865845, 140.269623403931, 513.231744616699, 309.207254891968,
    231.969481045532, 69.8056511077881, 112.592157028198, 133.99590116272,
    122.84378286438, -121.360155963135, 142.463153100586, 160.598022171021,
    85.8361654632568, -115.151129296875, 168.646501400757, 486.44165673523,
    377.756587875366, 53.8957852569581, -169.407985803223, 44.9826584106446,
    443.977754425049, 275.398571392822, -310.315822311401, -43.3600234497069,
    311.147247875977, 113.807384719849, -177.500178259277, 65.3731299591064,
    54.9020637268067, 570.606228424072, 550.219253933716, 203.623202307129,
    269.51739392395, -533.095908920288, 633.813805151368, 41.3413466583252,
    61.4782814941406, 27.2753046936034, -352.979056933594, 42.8206022369385,
    0.023822993829017, 0.025214359832736, 2.19169732505354, 37.5264082272386,
    15.364486108151, 15.509106911437, 9.61508703457049, 6.1518878216691,
    4.3628672426463, 3.3802552873144, 3.22763945154522, 1.70404365345812,
    6.75124076685407, 2.52192468892729, 5.32909357673836, 3.45829401531024,
    7.98970961672677, 19.9336693547378, 0.0293138809489544, 0.0329917711988839,
    0.0250384655566095, 0.0281001330035308, 0.0413807803554859,
    0.0320538360687485, 23.163075179387, 23.0248334748196, 5.50651705913325,
    9.27445399002468, 4.97179872424223, 2.0332891321427, 17.9817420667401,
    11.6915023984785, 1.21323115679805, 0.66298113784911, 0.758228735297754,
    0.346999668468278, 5.21196213017555, 49.3425691216675, 55.6319957709478,
    5.98373706386212, 1.20466151640824, 0.569844300791712, 0.592282774780343,
    0.364085619681276, 0.46390836075933, 0.592595071189114, 0.546624086593939,
    0.451476944683113, 0.584429675377772, 0.641752996871491, 2.4911949119413,
    0.528820306502942, 0.567834564975332, 4.88674382005111, 4.21575255159655,
    0.342986475466461, 0.597423022749927, 0.99662120299691, 5.78665756329068,
    2.74039503855587, 4.87164984176484, 1.35008109397666, 2.02326196743239,
    3.07320566450009, 1.78386143039334, 0.587567537122634, 1.0253510967396,
    0.579926783486469, 0.563366142039249, 2.18828874851807, 0.568788425490468,
    1.50505225904744, 0.497841164221427, 3.1035317403357, 4.35482051545463,
    4.53647101259057, 11.8846568027383, 1.66748060418872, 0.944403626857751,
    10.9183271600232, 10.8574785862895, 34.9902278482295, 23.930857949513,
    16.2234808928114, 3.37041675376762, 1.81767101188337, 3.88488155994652,
    1.91853479301712, 14.3992573317812, 39.5536245711944, 42.8593773147572,
    38.5487130634045, 38.8481439487416, 39.5447591753502, 38.0928464661996,
    35.4185354327967, 35.7541860516039, 32.8947672895871, 32.9898801259725,
    34.7393628296198, 34.4184868738917, 35.83742016615, 35.8010227010442,
    37.0158071277016, 35.1949936065939, 32.8535772823041, 33.4055116103626,
    32.8036973330392, 30.7208660565316, 31.1471238229436, 35.4102532658925,
    35.1837687094704, 39.2670197099304, 41.8351758162839, 40.1183898180856,
    36.3029394906325, 36.1961776319097, 41.2505176565543, 46.1172924990858,
    44.1624472537161, 41.7188537006873, 41.6223344578442, 38.7037348538098,
    40.9761115210062, 35.3348585118595, 38.6973063951532, 37.8845986725255,
    33.5468787081693, 35.5856462465408, 34.5168101225047, 37.5019839698826,
    33.4122316608416, 35.2894140982971, 38.1114786701591, 37.5594259670455,
    38.2652525483451, 40.2069488877486, 39.2440739772946, 39.409479205415,
    40.2402617843081, 41.1395358671476, 44.8776630216142, 41.4139692070141,
    39.0821712890854, 42.4900522367758, 39.6591700192056, 35.1116779843278,
    35.0522939188049, 37.864766216241, 38.1343736157423, 37.1858796166918,
    38.8131564957105, 35.2741423707838, 39.3939792375241, 37.8997578615876,
    38.4075339736066, 43.0367101785036, 36.2679418096117, 38.2498692137018,
    40.4662139666954, 36.2071212022431, 39.6788744943134, 37.2763955547187,
    37.6268625154569, 36.938169836219, 42.8989065074945, 39.1688267735901,
    41.4746585764728, 37.0780032991935, 37.3420213753244, 38.3742548295513,
    35.5753961044959, 38.2553481588851, 38.1194073751329, 37.5827760360015,
    37.2871511654341, 40.480395035647, 37.2062796404646, 40.2328318571221,
    37.2174175393717, 37.4562457440803, 38.6579017130362, 38.8312383516646,
    38.4138239749506, 37.9559840471282, 38.9627926002009, 37.7533876962371,
    38.7575158368563, 79.7419286483789, 75.2170638672125, 86.4778485619392,
    96.5800122108511, 87.1955294213962, 97.4236897024906, 92.0500612718797,
    94.7898696971784, 93.932276582702, 90.9616071082013, 91.8001587964608,
    90.8751804765465, 84.7559288122112, 92.8115222442691, 98.857488642662,
    96.5326549488763, 92.0468217869854, 90.9089203023527, 95.516373768878,
    88.6575256892371, 93.9187010168014, 91.9875851581135, 94.3459614967543,
    92.2897103061876, 92.4423320862505, 96.9779969347479, 95.8601271433708,
    96.2933656243666, 97.0628867868788, 93.6528420434592, 96.0180970921788,
    98.4787590138852, 101.641094541226, 93.5078861812489, 98.1703143662949,
    103.430612294944, 105.283446414585, 110.734886731358, 103.867242803756,
    106.644336590468, 101.597001893097, 108.656136975888, 93.026234965923,
    104.96257892017, 114.507353792363, 113.874808054894, 113.762365070777,
    113.574262628089, 105.881038659388, 105.833740223785, 110.286003082992,
    114.783353500622, 107.237564965082, 117.757336354112, 148.498810522795,
    133.418767027364, 156.185513071828, 0.0542044006347656, 0.0778095428466797,
    14.7540881469727, 141.873898809814, 35.4042162597656, 50.081369128418,
    7.82379324645996, 17.8349963378906, 5.84096130065918, 13.5170038421631,
    8.19710419921875, 4.88451590881348, 22.733500479126, 4.25417118530273,
    25.1368536621094, 8.85192832946777, 29.3840307312012, 20.1675340942383,
    0.100540420532227, 0.0778095428466797, 0.0603242523193359,
    0.0375933746337891, 0.153870556640625, 0.107534536743164, 111.289502883911,
    78.6934242828369, 10.9562830444336, 15.7638636749268, 11.493081463623,
    3.95692124633789, 37.7201429901123, 41.1498827270508, 1.32013943481445,
    1.22659313049316, 3.72349261779785, 0.597122671508789, 5.02964382019043,
    170.532289984131, 61.917162286377, 24.4199567504883, 3.42012282714844,
    1.06397992858887, 1.57717320556641, 1.22659313049316, 0.994038766479492,
    1.66984524536133, 0.806946157836914, 0.783341015625, 1.13392109069824,
    2.34652598876953, 6.77992140197754, 0.690668975830078, 0.993164501953125,
    6.7091059753418, 9.39397233581543, 0.783341015625, 2.47416860961914,
    3.19805963745117, 21.9807587219238, 13.7679177612305, 10.1974214355469,
    3.22166477966309, 8.51795928039551, 3.52416030578613, 5.46153049621582,
    1.63487466430664, 4.10816900939941, 1.63487466430664, 1.05174022521973,
    7.95755571899414, 0.67755500793457, 3.9674124206543, 1.37784089355469,
    4.90112693481445, 17.8979433837891, 12.9513546936035, 63.0711914611816,
    7.53790874633789, 3.92107640075684, 32.9011969207764, 25.7610785339355,
    106.260733328247, 96.5773794342041, 50.3078036407471, 9.94126192932129,
    2.61492519836426, 8.37720269165039, 3.99101756286621, 60.1310398590088,
    89.248419909668, 97.1351602020264, 79.3552425292969, 77.8847295959473,
    80.6613937316895, 86.9849490509033, 78.4914691772461, 70.7218803314209,
    66.498308404541, 67.2912663299561, 79.6813431976318, 76.9046790618897,
    74.5485361633301, 77.4414774810791, 84.6043267456055, 88.8051677947998,
    86.5880329559326, 75.0853345825195, 96.9489418579102, 83.8113688201904,
    75.691199899292, 85.6551927062988, 84.81502449646, 91.1377055511475,
    92.3048486938477, 82.3417301513672, 79.4479145690918, 73.5212753448486,
    92.1649663696289, 98.3250342224121, 92.4683361602783, 90.2048653015137,
    95.3149414581299, 86.9386130310059, 82.7386462463379, 85.9585624969483,
    88.0113356048584, 78.4914691772461, 68.9716027496338, 67.80533387146,
    75.7847462036133, 85.3509486511231, 59.6379546661377, 82.5751587799072,
    86.5416969360352, 82.6914359619141, 79.7512843597412, 83.4616630096436,
    73.5684856292725, 75.2715529266357, 92.5854876068115, 98.8382274993897,
    94.8253533233643, 90.2048653015137, 88.3383105377197, 88.0812767669678,
    91.0677643890381, 79.7285534820557, 67.2449303100586, 85.4217640777588,
    85.5616464019775, 77.8147884338379, 78.3052508331299, 81.9212089141846,
    84.8613605163574, 79.0754778808594, 93.8680336669922, 94.6618658569336,
    72.1451829803467, 83.7414276580811, 81.3380744750977, 75.2715529266357,
    77.3479311767578, 70.8154266357422, 71.0479809997559, 81.4552259216309,
    83.0647469146729, 78.6549566436768, 93.8216976470947, 80.8712172180176,
    68.2249808441162, 64.6090227630615, 71.4212919525147, 69.881712121582,
    88.3846465576172, 69.5783423309326, 79.8684358062744, 85.8178059082031,
    79.8684358062744, 83.998461428833, 73.4049981628418, 74.7120236297607,
    73.3114518585205, 81.9448140563965, 80.2181416168213, 84.7678142120361,
    79.3779734069824, 76.0181748321533, 88.8278986724854, 224.616042114258,
    220.533226776123, 233.390160900879, 270.862887030029, 196.312602337647,
    246.783019180298, 241.88276651001, 279.076602255249, 243.376884585571,
    214.513041247559, 234.532824636841, 276.813131396484, 229.049437527466,
    249.769506802368, 236.166825036621, 247.506035943604, 246.783019180298,
    240.529405023193, 247.482430801392, 244.892859274292, 244.449607159424,
    235.932522143555, 218.59673085022, 224.8031347229, 222.936579959106,
    243.352405178833, 227.556193716431, 259.966053973389, 282.599014031982,
    225.945798458862, 231.430059832764, 250.959380822754, 323.899270257568,
    229.913210879517, 236.959782962036, 244.846523254395, 241.929976794434,
    271.02637449646, 214.886352200317, 207.139494232178, 226.716899771118,
    316.01252996521, 201.119308703613, 239.05976635437, 260.712675878906,
    233.762597589111, 275.109189834595, 265.730079995728, 263.722768643189,
    289.249544284058, 225.012958209229, 226.693294628906, 267.082567218018,
    256.093062121582, 330.759623995972, 336.429229449463, 331.109329806519,
    -0.0620727813720703, -0.0620727813720703, -4.12215724182129,
    -74.3081134185791, -37.7210172546387, -38.0943282073975, -45.8184552978516,
    -16.7683936157227, -19.988309866333, -9.86170385742187, -7.69177930297852,
    -4.3083755859375, -22.0646881164551, -9.13781282958984, -15.1816035003662,
    -14.1552169464111, -15.5313093109131, -76.0575167358398, -0.0384676391601563,
    -0.0620727813720703, -0.0559529296875, -0.0568271942138672,
    -0.103163214111328, -0.0568271942138672, -1.01327258605957,
    -29.5719976043701, -15.7594923522949, -30.5284429962158, -16.4816348510742,
    -5.00254161987305, -40.655049005127, -33.2115608276367, -4.02161682128906,
    -1.73628934936523, -1.52646586303711, -0.826179977416992, -27.5655605163574,
    -85.0117340148926, -239.641152264404, -9.76203770141602, -4.44301232299805,
    -1.33937325439453, -1.10681889038086, -0.686297653198242, -0.989667443847656,
    -1.31664237670898, -2.2966929107666, -0.732633673095703, -1.5028607208252,
    -1.19949093017578, -8.94547463378906, -1.87704593811035, -1.6200121673584,
    -17.4389545074463, -8.52495339660645, -1.05960860595703, -1.37521809997559,
    -3.26450374145508, -12.9487319000244, -10.9414205474854, -13.2975634460449,
    -6.08837816162109, -3.70775585632324, -10.0785214599609, -5.85407526855469,
    -1.56231070861816, -1.23533577575684, -1.49149528198242, -1.84207535705566,
    -4.26815941772461, -3.21816772155762, -5.5979157623291, -1.32888208007813,
    -19.807337109375, -9.96049574890137, -10.9178154052734, -18.3840344604492,
    -3.31171402587891, -1.7957393371582, -30.3544643554688, -36.8179019989014,
    -139.180289804077, -58.5643578277588, -31.9176493286133, -11.664437310791,
    -8.81783201293945, -15.2576645141602, -5.99483185729981, -29.0011028686523,
    -81.6904030792236, -91.0931180603027, -79.8929152130127, -78.1662427734375,
    -73.8963348266602, -60.4798714050293, -71.9834440429688, -58.286341708374,
    -59.4534848510742, -62.1838129669189, -55.5804929992676, -58.6832578033447,
    -69.0432924407959, -57.3762323364258, -59.3136025268555, -47.1106182678223,
    -54.576837322998, -59.7804597839356, -48.4167694702148, -40.3202056915283,
    -50.3768705383301, -64.1666449127197, -57.2136191345215, -68.7399226501465,
    -73.5238981384277, -68.5300991638184, -66.5699980957031, -50.4704168426514,
    -69.5564857177734, -76.3232931518555, -74.7137721588135, -69.3938725158691,
    -65.0994851623535, -77.7238649230957, -80.8738400115967, -62.5107878997803,
    -74.3867972259522, -71.5401919281006, -60.7605103179932, -49.7701309570313,
    -61.5062579589844, -60.7832411956787, -72.8699482727051, -64.1902500549316,
    -57.5169889251709, -67.4337714477539, -65.7070990081787, -54.60044246521,
    -69.8371246307373, -67.6199897918701, -49.3732148620606, -65.5208806640625,
    -68.5073682861328, -62.7660731414795, -54.7403247894287, -80.4533187744141,
    -78.7039154571533, -65.3565189331055, -61.6470145477295, -69.9534018127442,
    -84.0701511199951, -67.8070824005127, -64.8896616760254, -52.5231899505615,
    -70.4202590698242, -54.2970726745605, -61.8568380340576, -82.4597558624268,
    -79.3805962005615, -77.2797385437012, -75.3432426177979, -77.3496797058105,
    -99.6801442382812, -96.3666816833496, -66.1966871429443, -66.0804099609375,
    -85.3771765869141, -86.6361175048828, -86.519840322876, -73.686511340332,
    -72.9862254547119, -97.6264968658447, -67.643594934082, -82.3434786804199,
    -76.2297468475342, -62.486308493042, -66.8969730285645, -73.4766878540039,
    -74.6665618743896, -63.3273509674072, -76.7665452667236, -77.6067134765625,
    -80.1036129638672, -75.6702175506592, -68.7399226501465, -71.6564691101074,
    -70.6300825561524, -47.2741057342529, -66.5472672180176, -119.851175390625,
    -115.464990261841, -165.187910934448, -151.211043951416, -144.747606307983,
    -177.087525402832, -145.120917260742, -184.741711331177, -178.791466964722,
    -149.927623626709, -166.237902630615, -154.897817459106, -138.845446490479,
    -188.824526669312, -161.267708798218, -181.077668701172, -186.351232324219,
    -166.120751184082, -159.144120263672, -147.290841815186, -152.657951742554,
    -146.754917660522, -149.881287606812, -159.401154034424, -157.30117064209,
    -163.928095751953, -173.400751895142, -148.084674005127, -170.460600292969,
    -146.381606707764, -154.687993972778, -197.130913934326, -179.420937423706,
    -165.864591677856, -184.064156323242, -183.108585195923, -173.471567321777,
    -197.644107211304, -188.917198709106, -196.034586218262, -192.394148730469,
    -162.621070285034, -166.611213583374, -172.887558618164, -208.331116781616,
    -179.724307214355, -212.904394519043, -216.054369607544, -173.961155456543,
    -198.927527536011, -195.357031210327, -206.371015713501, -202.0538974823,
    -206.604444342041, -214.583856674194, -232.924177908325, -221.771185345459,
    0.410030062866211, 0.421395501708984, 13.8400445846558, 563.359449765015,
    33.6185309646606, 12.7393455459595, -93.0300511184693, -8.25873984832764,
    -65.091179649353, 4.21614067840577, 0.302495526123048, -23.6667778610229,
    22.1223895751953, -35.2949331939697, 19.050224029541, 0.381179333496095,
    35.9357690917969, -376.63665501709, 1.01545824737549, 0.282387442016602,
    -0.17878709564209, -0.480845489501953, -0.0987918914794922,
    0.114091520690918, 309.479588291931, 336.861553257751, -15.8849493118286,
    -102.701602441406, 39.7987069015503, -29.943122895813, -3.00659570617675,
    52.3216719772339, -30.2473669509888, -7.05138053741455, 2.79852074890137,
    -2.81250898132324, -65.2231935928345, 920.133689007568, -590.602843803406,
    -0.685423388671879, -24.5257427581787, -3.57049632568359, 4.06795284118652,
    1.83814116668701, -0.35276573638916, 0.370688159179688, -12.7292915039063,
    -1.71924119110107, -1.55487946014404, 1.21959901428223, 2.26609365234375,
    -20.4800836624146, -6.10761198120117, -68.6420050231934, -9.67461124877929,
    1.70219303283691, -1.72011545562744, 3.61114962615967, 35.7093345794678,
    1.35685854492187, -3.79518230895996, -25.5320212280273, 18.5326594299316,
    -85.8672018539429, -11.9472618850708, -16.1900676315308, 10.2481287780762,
    -8.1157975982666, -5.0943393951416, 21.8247025039673, -12.8285205276489,
    -5.89735136260987, -3.25619822845459, -4.9050611251831, 5.1555379119873,
    43.7962814483643, 131.995146794128, -1.63837172241211, 5.1502923248291,
    74.1284520584106, 16.0545566299438, -245.215900016785, -83.0520700790405,
    102.326542959595, -5.0331408782959, -5.28973751678467, 2.1874098449707,
    10.5169651199341, -75.3327514434815, 293.693430871582, 419.047227191162,
    -73.2839125259399, -30.562102180481, 182.712980497742, 232.317438327026,
    -211.587315010071, 210.669337257385, -119.435025476074, 268.449916937256,
    236.731162788391, 420.957495181274, -172.387916441345, 267.928418147278,
    79.4706454467774, 341.192659721374, 221.042048730469, 242.423936251831,
    381.524230851746, 86.3790837341308, 226.325666395569, -236.441781230164,
    99.5813523468017, 777.705506488037, 289.040157929993, -0.322166477966277,
    -89.6654440887451, 277.645431225586, -236.006397496033, 385.779713433838,
    -192.385843217468, 317.980499414063, 381.872625265503, -51.46707840271,
    200.353015846252, -191.684245935059, 230.052218939209, -10.614445614624,
    -29.3009756011963, 204.286331950378, -60.2735449768066, 409.592056338501,
    32.9497186019898, 506.357402645874, 581.266135794067, 89.5539753616333,
    400.567897897339, 199.587160121155, -201.135482597351, 84.7844252380371,
    194.33720164032, -115.062828579712, 231.736052416992, 180.586769169617,
    -216.419812179565, 1.41456000366212, 337.292565669251, -333.452358737183,
    5.74566646728511, 190.329573051453, -123.611824250793, -37.8088808395386,
    265.542987387085, 4.47361158142095, 381.038576907349, 148.838290026855,
    -214.474573608398, 109.14318347168, 251.593659736633, -168.994458682251,
    526.183536442566, 190.834897947693, -172.314041088867, 95.0727701843262,
    196.925024638367, -234.690192251587, 605.992522261047, 211.891559065247,
    95.210466847229, -154.320365739441, 40.2039285095215, 214.869304042053,
    153.862688259888, -50.6199160766602, -57.1602889984131, 221.103684379578,
    171.900076835632, -94.0761086242676, -92.63881774292, 180.286896437073,
    179.703324865723, -84.7809281799316, 1.37696662902835, 269.113920845032,
    198.741309191895, -61.6859193191529, 64.9241951248169, 256.457630429077,
    -286.601834165955, 6.03504802551274, 755.754035627747, 183.450859757996,
    131.637135470581, 102.157372773743, 67.1837317932129, 193.882584086609,
    165.866777339172, 182.237817727661, -203.436983963013, -79.3705421585083,
    330.345659742737, 392.353308407593, 133.086666055298, 26.8386095626831,
    153.269499778748, -40.0343211914063, -2.72377113189694, 193.469931230164,
    535.107154463196, 407.171654997253, 286.044053398132, 81.6711692596435,
    100.704782263184, 123.60177020874, 131.792754556274, -139.465737171936,
    142.592107118225, 216.317960362244, 128.862656996155, -116.783381167603,
    156.899008959961, 473.807660064697, 345.849429721069, 41.0869356811524,
    -61.1762231002807, 183.815865197754, 331.75759695282, 232.186298648071,
    -256.217207684326, -31.9036610961914, 324.856152781677, 92.7848199188232,
    -138.008338206482, 196.69203314209, 179.490004321289, 463.799516899109,
    407.525294998169, 225.545385305786, 316.032638049317, -416.114943969726,
    494.574065359497, 84.9868174758911, 72.8660140823364, 109.420762458801,
    -236.557621279907, 144.039889173889, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 3.0, 5.0, 2.0, 3.0, 4.0, 3.0, 2.0, 4.0, 2.0,
    3.0, 5.0, 3.0, 2.0, 3.0, 3.0, 3.0, 3.0, 4.0, 4.0, 2.0, 3.0, 2.0, 5.0, 7.0,
    5.0, 4.0, 3.0, 3.0, 3.0, 5.0, 3.0, 5.0, 5.0, 3.0, 3.0, 2.0, 4.0, 3.0, 2.0,
    5.0, 4.0, 7.0, 3.0, 4.0, 6.0, 4.0, 5.0, 3.0, 3.0, 2.0, 4.0, 4.0, 6.0, 5.0,
    3.0, 3.0, 5.0, 2.0, 3.0, 4.0, 5.0, 3.0, 3.0, 4.0, 6.0, 4.0, 2.0, 4.0, 4.0,
    3.0, 5.0, 3.0, 2.0, 4.0, 4.0, 3.0, 3.0, 3.0, 5.0, 2.0, 2.0, 4.0, 3.0, 3.0,
    3.0, 5.0, 3.0, 3.0, 2.0, 4.0, 3.0, 4.0, 3.0, 4.0, 4.0, 2.0, 3.0, 5.0, 4.0,
    9.0, 7.0, 6.0, 6.0, 4.0, 4.0, 4.0, 5.0, 6.0, 5.0, 4.0, 6.0, 5.0, 6.0, 7.0,
    6.0, 5.0, 5.0, 6.0, 7.0, 5.0, 4.0, 5.0, 4.0, 4.0, 6.0, 4.0, 5.0, 6.0, 4.0,
    4.0, 6.0, 5.0, 6.0, 6.0, 5.0, 4.0, 5.0, 6.0, 5.0, 5.0, 4.0, 6.0, 7.0, 6.0,
    4.0, 6.0, 7.0, 5.0, 6.0, 5.0, 7.0, 7.0, 6.0, 7.0, 5.0, 7.0, 0.0, 0.0, 0.0,
    10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 11.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 36.0, 38.0, 18.0,
    19.0, 36.0, 36.0, 20.0, 22.0, 20.0, 23.0, 29.0, 41.0, 19.0, 40.0, 23.0, 38.0,
    39.0, 23.0, 22.0, 17.0, 39.0, 22.0, 40.0, 27.0, 42.0, 37.0, 19.0, 38.0, 22.0,
    42.0, 24.0, 40.0, 39.0, 19.0, 38.0, 19.0, 38.0, 38.0, 19.0, 41.0, 22.0, 43.0,
    21.0, 41.0, 37.0, 40.0, 39.0, 37.0, 19.0, 19.0, 38.0, 24.0, 40.0, 39.0, 20.0,
    34.0, 41.0, 20.0, 18.0, 40.0, 37.0, 20.0, 40.0, 24.0, 38.0, 38.0, 18.0, 36.0,
    37.0, 21.0, 37.0, 36.0, 18.0, 38.0, 42.0, 19.0, 33.0, 36.0, 39.0, 18.0, 18.0,
    37.0, 37.0, 34.0, 22.0, 38.0, 34.0, 36.0, 18.0, 37.0, 36.0, 36.0, 32.0, 37.0,
    37.0, 19.0, 19.0, 38.0, 36.0, 39.0, 40.0, 38.0, 39.0, 37.0, 38.0, 39.0, 38.0,
    40.0, 26.0, 25.0, 35.0, 41.0, 39.0, 39.0, 36.0, 27.0, 26.0, 39.0, 40.0, 40.0,
    38.0, 39.0, 38.0, 38.0, 38.0, 39.0, 37.0, 38.0, 38.0, 37.0, 39.0, 38.0, 40.0,
    39.0, 36.0, 36.0, 46.0, 38.0, 36.0, 37.0, 37.0, 38.0, 38.0, 38.0, 36.0, 36.0,
    37.0, 36.0, 37.0, 36.0, 45.0, 36.0, 38.0, 44.0, 35.0, 41.0, 0.0, 0.0, 0.0,
    5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 2.75, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 12.0, 7.6, 9.0,
    6.33333333333333, 9.0, 12.0, 10.0, 5.5, 10.0, 7.66666666666667, 5.8,
    13.6666666666667, 9.5, 13.3333333333333, 7.66666666666667, 12.6666666666667,
    13.0, 5.75, 5.5, 8.5, 13.0, 11.0, 8.0, 3.85714285714286, 8.4, 9.25,
    6.33333333333333, 12.6666666666667, 7.33333333333333, 8.4, 8.0, 8.0, 7.8,
    6.33333333333333, 12.6666666666667, 9.5, 9.5, 12.6666666666667, 9.5, 8.2,
    5.5, 6.14285714285714, 7.0, 10.25, 6.16666666666667, 10.0, 7.8,
    12.3333333333333, 6.33333333333333, 9.5, 9.5, 6.0, 6.66666666666667, 7.8,
    6.66666666666667, 11.3333333333333, 8.2, 10.0, 6.0, 10.0, 7.4,
    6.66666666666667, 13.3333333333333, 6.0, 6.33333333333333, 9.5, 9.0, 9.0,
    9.25, 7.0, 7.4, 12.0, 9.0, 9.5, 10.5, 6.33333333333333, 11.0, 12.0, 7.8, 9.0,
    9.0, 9.25, 12.3333333333333, 11.3333333333333, 7.33333333333333, 7.6,
    11.3333333333333, 12.0, 9.0, 9.25, 12.0, 9.0, 10.6666666666667, 9.25, 9.25,
    9.5, 6.33333333333333, 7.6, 9.0, 4.33333333333333, 5.71428571428571,
    6.33333333333333, 6.5, 9.25, 9.5, 9.75, 7.6, 6.66666666666667, 5.2, 6.25,
    5.83333333333333, 8.2, 6.5, 5.57142857142857, 6.0, 5.4, 5.2, 6.5,
    5.71428571428571, 8.0, 9.5, 7.8, 9.5, 9.5, 6.33333333333333, 9.75, 7.4,
    6.33333333333333, 9.5, 9.25, 6.5, 7.6, 6.66666666666667, 6.5, 7.2, 9.0, 9.2,
    6.33333333333333, 7.2, 7.4, 9.25, 6.33333333333333, 5.42857142857143,
    6.33333333333333, 9.0, 6.0, 5.28571428571429, 7.2, 6.16666666666667, 7.2,
    6.42857142857143, 5.14285714285714, 6.33333333333333, 6.28571428571429, 7.0,
    5.85714285714286, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 1.0, 2.0, 4.0,
    2.0, 4.0, 4.0, 3.0, 3.0, 2.0, 4.0, 3.0, 2.0, 2.0, 4.0, 4.0, 2.0, 3.0, 4.0,
    3.0, 4.0, 3.0, 2.0, 3.0, 3.0, 3.0, 2.0, 3.0, 2.0, 4.0, 3.0, 3.0, 3.0, 2.0,
    3.0, 3.0, 4.0, 4.0, 2.0, 4.0, 4.0, 4.0, 2.0, 5.0, 4.0, 4.0, 4.0, 4.0, 3.0,
    3.0, 1.0, 5.0, 4.0, 5.0, 10.0, 5.0, 10.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 2.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 6.0, 4.0, 3.0, 6.0, 2.0, 2.0, 4.0, 4.0, 6.0, 5.0, 4.0, 7.0, 8.0,
    4.0, 1.0, 6.0, 2.0, 0.0, 4.0, 4.0, 6.0, 9.0, 5.0, 8.0, 7.0, 8.0, 6.0, 8.0,
    4.0, 8.0, 8.0, 5.0, 8.0, 6.0, 8.0, 5.0, 9.0, 11.0, 7.0, 5.0, 8.0, 8.0, 6.0,
    9.0, 9.0, 8.0, 8.0, 8.0, 9.0, 11.0, 9.0, 7.0, 7.0, 5.0, 7.0, 5.0, 50.0, 50.0,
    50.0, 48.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 49.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 44.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 49.0, 49.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 49.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 49.0, 50.0, 50.0, 50.0, 50.0, 50.0, 49.0, 50.0, 50.0, 50.0, 49.0, 49.0,
    49.0, 50.0, 48.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 49.0, 49.0, 48.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 49.0, 49.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 49.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 46.0, 43.0, 44.0, 43.0, 42.0, 44.0, 44.0, 43.0, 43.0,
    42.0, 41.0, 43.0, 41.0, 40.0, 42.0, 45.0, 42.0, 45.0, 46.0, 43.0, 42.0, 41.0,
    39.0, 42.0, 39.0, 40.0, 40.0, 41.0, 40.0, 42.0, 39.0, 39.0, 42.0, 40.0, 41.0,
    39.0, 41.0, 37.0, 37.0, 39.0, 41.0, 38.0, 40.0, 39.0, 37.0, 37.0, 38.0, 38.0,
    39.0, 38.0, 38.0, 36.0, 39.0, 38.0, 35.0, 38.0, 35.0, 0.0236051422119141,
    0.0, 98.4194547912598, 99.7028751159668, 1.18987402038574, 2.31068114318848,
    10.6389250213623, 17.1749266204834, 0.302495526123048, 3.40613459472656,
    1.19162254943848, 1.44690779113769, 9.82323621826172, 13.3465222595215,
    1.12080712280273, 4.66594977722168, 5.99570612182617, 56.1164171539307,
    0.0236051422119141, 0.0437132263183594, 0.0, 0.0, 0.326974932861328,
    0.280638912963867, 25.8065402893066, 42.0460038665772, 1.91289078369141,
    3.14997508850098, 4.76037034606934, 0.349705810546874, 29.7284909545898,
    11.0594462585449, 0.535924154663086, 0.769352783203125, 0.232554364013672,
    0.092672039794922, 13.1603039154053, 8.81870627746582, 149.706434701538,
    2.61405093383789, 0.792083660888672, 0.023605142211914, 0.0708154266357421,
    0.0, 0.0463360198974609, 0.39779035949707, 0.0463360198974609,
    0.140756588745117, 0.162613201904297, 0.886504229736328, 0.18621834411621,
    0.140756588745117, 0.187092608642578, 2.63590754699707, 1.2842945892334,
    0.00087426452636713, 1.09720198059082, 0.653075601196289, 2.37887377624512,
    1.84294962158203, 3.87386611633301, 0.467731521606444, 0.396041830444336,
    3.70950438537598, 6.25273989257813, 0.932840249633789, 0.093546304321289,
    0.304244055175781, 0.093546304321289, 2.56684064941406, 0.0463360198974612,
    2.49689948730469, 0.770227047729492, 8.58527764892578, 0.887378494262695,
    4.64321889953613, 7.06930296020508, 3.07915966186524, 2.52137889404297,
    1.75027758178711, 0.560403561401367, 6.41622735900879, 9.56620244750977,
    32.4561962768555, 8.4008078338623, 2.70672297363281, 3.6850249786377,
    3.0800339263916, 10.8732279144287, 2.03091649475098, 5.94937010192872,
    34.3708355895996, 11.1765977050781, 11.3855469268799, 7.02296694030763,
    17.5464890441895, 10.1484626220703, 9.7305641784668, 22.6583137298584,
    21.0706493499756, 12.5054797851562, 17.6863713684082, 8.58527764892578,
    1.68033641967774, 9.26283265686035, 3.40700885925293, 1.16714314270018,
    8.70505188903809, 5.43617682495118, 4.94658869018554, 12.0158916503906,
    9.75329505615235, 5.2490842163086, 9.98672368469238, 19.8790268005371,
    13.790648638916, 0.464234463500972, 19.902631942749, 14.6544219909668,
    4.19909252014161, 30.402548904419, 43.5173910644531, 22.888245300293,
    6.53337880554199, 7.49069846191406, 2.42608406066894, 3.92020213623047,
    11.1075308074951, 2.70497444458008, 8.28278212280274, 15.1894718811035,
    3.98839476928711, 21.8181455200195, 12.9260010223389, 8.60975705566408,
    6.46343764343262, 12.1557739746094, 25.9228174713135, 14.2321522247314,
    8.8912702331543, 2.75393325805663, 23.1706327423096, 11.0830514007568,
    11.20020284729, 9.28818632812499, 7.02384120483399, 13.5554714813232,
    15.65720340271, 0.0699411621093731, 26.3888004638672, 6.62605084533692,
    0.419646972656253, 8.86679082641602, 1.75027758178712, 22.5866240386963,
    2.77491560668945, 5.10920189208984, 9.14655547485351, 3.64043748779297,
    19.109674017334, 2.82387442016601, 5.83396718444824, 19.1332791595459,
    13.0431524688721, 16.450161328125, 0.676680743408213, 10.7097404479981,
    32.3635242370605, 26.3424644439697, 25.0354389770508, 15.8206908691406,
    1.70306729736328, 0.163487466430674, 4.59688287963867, 7.46621905517578,
    6.41622735900879, 7.83953000793458, 9.14568121032715, 1.98370621032716,
    8.47162326049805, 16.3557407592773, 0.326974932861333, 5.13193276977539,
    12.1094379547119, 7.6760425415039, 16.1004555175781, 12.3437408477783,
    24.336027355957, 42.9788441162109, 22.8200526672363, 40.4364828735351,
    64.1439140350342, 25.5022962341309, 39.3628860351563, 27.2761789581299,
    17.4057324554444, 42.9797183807373, 69.3929982513428, 41.8597855224609,
    12.6226312316894, 34.8831546020508, 3.73310952758789, 13.8588412719727,
    47.4367189361572, 34.7441465423584, 56.443392086792, 78.6086206237793,
    69.8353761016846, 67.2702839813233, 85.3763023223877, 29.2817417816162,
    71.8654183319092, 58.6823835388184, 65.145821182251, 13.8605898010254,
    57.7040815338135, 53.6894588287354, 71.8418131896973, 46.0833574493408,
    48.8818781982422, 154.513141067505, 12.9277495513916, 29.7267424255371,
    0.561277825927704, 111.626094726562, 75.6684690216064, 64.2593169525147,
    41.0414739257813, 11.8760093261719, 18.0360771789551, 42.1168192932129,
    70.6056031494141, 48.4167694702149, 8.09656377868652, 117.995986065674,
    26.7411290679932, 111.227430102539, 22.3059851257324, 9.21737090148926,
    41.6508363006592, 37.0285997497559, 46.9934668212891, 24.2433553161621,
    59.9658038635254, 19.2032203216552, -0.000548956795625908, 0.0,
    98.4194547912598, -49.8514375579834, 1.18987402038574, -0.770227047729492,
    3.5463083404541, 5.72497554016113, 0.151247763061524, -1.70306729736328,
    0.397207516479492, -0.361726947784424, 9.82323621826172, 13.3465222595215,
    0.0622670623779297, 4.66594977722168, 5.99570612182617, 28.0582085769653,
    0.0236051422119141, -0.00230069612201891, 0.0, 0.0, 0.040871866607666,
    0.280638912963867, 12.9032701446533, 42.0460038665772, 0.159407565307617,
    1.57498754425049, 1.19009258651733, -0.0388562011718749, 3.30316566162109,
    2.21188925170898, 0.107184830932617, -0.769352783203125, 0.00930217456054688,
    -0.00545129645852482, -0.626681138828823, -4.40935313873291,
    49.9021449005127, -0.153767701990464, 0.158416732177734,
    -0.00393419036865234, 0.0354077133178711, 0.0, -0.0115840049743652,
    0.0137169089481748, 0.00421236544522372, 0.140756588745117,
    0.0147830183549361, 0.221626057434082, 0.18621834411621, -0.140756588745117,
    0.187092608642578, -0.658976886749268, 0.0987918914794922,
    0.000218566131591783, -1.09720198059082, -0.21769186706543, 0.59471844406128,
    -1.84294962158203, -0.297989701256385, 0.233865760803222,
    -0.0282887021745954, -1.85475219268799, 0.694748876953125, 0.932840249633789,
    -0.0311821014404297, -0.101414685058594, -0.093546304321289,
    -2.56684064941406, -0.0463360198974612, -1.24844974365234, 0.770227047729492,
    -4.29263882446289, 0.887378494262695, -4.64321889953613, -3.53465148010254,
    -1.53957983093262, 0.840459631347656, -0.0673183685302735,
    -0.560403561401367, 6.41622735900879, 9.56620244750977, -16.2280981384277,
    -2.8002692779541, -0.112780123901367, 0.118871773504442, -0.513338987731934,
    10.8732279144287, 0.290130927821568, 0.424955007280623, 2.29138903930664,
    0.385399920864763, 0.517524860312722, -0.234098898010254, 17.5464890441895,
    2.02969252441406, 0.3355366958092, -3.23690196140834, -1.7558874458313,
    0.893248556082589, -0.680245052631085, 0.953919738769531, 0.0560112139892579,
    0.544872509227079, -0.121678887830462, -0.0376497787967801,
    -0.725420990753174, -0.543617682495118, -1.64886289672851, 12.0158916503906,
    0.886663186922941, -0.194410526529948, -0.369878654988607, 2.20878075561524,
    -1.97009266270229, 0.154744821166991, -0.663421064758301, 0.505324896240235,
    0.59987036002023, 10.1341829681397, -3.34749162034255, 0.738330493557838,
    -6.53337880554199, -0.468168653869629, 0.0836580710575497, 0.392020213623047,
    1.58679011535645, -0.0819689225630328, -0.345115921783447, 0.523774892451846,
    -0.30679959763747, 2.18181455200195, -0.497153885474572, 0.614982646833148,
    -0.215447921447754, -0.935059536508413, 12.9614087356567, 2.37202537078857,
    -0.404148646961559, 2.75393325805663, 4.63412654846191, 11.0830514007568,
    -11.20020284729, 0.371527453125, -0.585320100402832, 0.645498641967774,
    1.42338212751909, -0.0099915945870533, -26.3888004638672, 0.662605084533692,
    -0.209823486328126, -0.216263190888195, -0.0407041298090027,
    1.32862494345272, -0.0895134066674015, -0.189229699707031, 2.28663886871338,
    0.21414338163488, 2.38870925216675, -0.0806821262904575, 1.94465572814941,
    -19.1332791595459, 1.63039405860901, 0.365559140625, 0.0845850929260266,
    5.35487022399903, 5.39392070617676, -1.54955673199822, 12.5177194885254,
    5.27356362304688, -1.70306729736328, -0.00480845489501981, 0.153229429321289,
    1.2443698425293, 6.41622735900879, 0.290352963256836, -0.703513939255934,
    0.198370621032716, -0.847162326049805, 0.860828461014597, 0.163487466430666,
    0.513193276977539, 0.931495227285533, 0.639670211791992, -0.947085618681067,
    -12.3437408477783, 0.901334346516927, 2.86525627441406, 1.14100263336182,
    -2.12823594071238, -2.06915851725917, -0.579597641684792, -2.07173084395559,
    -3.89659699401855, 1.24326660396031, -2.26209044109144, 3.65226306586015,
    6.97663092041015, 2.10377187194824, 34.8831546020508, -0.196479448820415,
    -1.38588412719727, 2.37183594680786, 1.82863929170307, 9.40723201446533,
    13.1014367706299, -9.97648230024065, 11.2117139968872, 7.76148202930797,
    -14.6408708908081, -2.31823930102933, -3.08854650204307, -2.03580691194534,
    -0.729504726369757, 9.61734692230225, -2.82576099098607, 11.9736355316162,
    7.68055957489014, -2.57273043148643, -22.0733058667864, 2.15462492523193,
    4.95445707092285, -0.0295409382067212, -6.56624086626838, 1.84557241516113,
    -4.94302438096267, -1.32391851373488, -3.95866977539063, 0.949267219945005,
    -3.00834423522949, 70.6056031494141, 6.91668135288784, 0.261179476731823,
    -6.2103150560881, 1.21550586672696, 5.05579227738814, -1.06218976789202,
    -0.512076161193848, -6.94180605010986, 1.68311817044345, -7.83224447021485,
    -0.897902048746745, 1.87393137073517, -3.84064406433105, -0.0375933746337891,
    -0.179224227905273, 348.393539492798, -469.197663217163, 77.8008002014161,
    70.9859082183838, -95.5789693450928, 38.7640148345947, -13.5344891326904,
    35.7556705993652, 29.3822822021484, -16.3539922302246, -52.2469223602295,
    7.48457861022949, 27.6180163879394, -7.72587561950684, -5.32864228820801,
    -91.2426172943114, 0.222937454223633, 0.224685983276367, 0.278016119384766,
    0.349705810546876, 0.396916094970703, -0.142505117797852, 98.8600841125488,
    -109.602172348022, 120.155419445801, 48.2725158233643, -146.433188314819,
    -2.057144430542, 172.162793325806, -93.8208233825684, 31.4630317749023,
    12.3970709838867, 0.712525588989258, 4.7026688873291, -31.7471677459717,
    775.256691549683, -1067.93684983521, 33.5962372192383, 4.27777632751465,
    -1.31926517028809, 0.382927862548828, 1.75814596252441, 1.82896138916016,
    -0.388173449707031, 3.83277568359375, 1.85518932495117, 2.6245421081543,
    -11.4012836883545, -15.5995019439697, 1.25107253723145, 5.32951655273438,
    11.9765497467041, 9.89754870300293, 2.76355016784668, -1.37609236450195,
    2.65951268920898, 37.7428738677978, -1.91988489990235, -32.3005771911621,
    -31.8774331604004, -3.25138977355957, 2.48378551940918, -8.90001287841797,
    6.64703319396973, 0.78683807373047, 0.556032238769531, -2.22063189697266,
    7.57200506286621, 6.20290681457519, 1.29915708618164, 1.37346957092285,
    14.366788961792, -47.2872197021484, -7.90335131835938, -51.0185807006836,
    39.5875720184326, -13.3771215179443, 3.57486764831543, -24.9969713378906,
    -157.330895635986, 72.027157269287, 91.5792091369629, 27.8566906036377,
    63.7627347015381, 73.755578237915, 24.7180809539795, -62.3385577880859,
    205.757282015991, 1317.78766323853, -153.553198617554, 110.8348853302,
    109.013792321777, 117.183794320679, -341.378440933228, 309.138187994385,
    -303.999261108398, 346.453546508789, 1175.15839987793, 768.435679714966,
    -220.953748013306, 289.544171429443, 771.374082788086, 639.282329763794,
    -27.8505707519531, 1242.0719839325, 1039.95776219788, -6.9696368041992,
    12.2344577819824, -798.499888247681, -352.74650256958, 2964.77526255798,
    308.93710715332, -219.62923725586, 25.1761955657959, 141.261913641357,
    -408.973951318359, 412.042619805908, -429.419501531982, 615.664073583984,
    727.627634417725, -32.8496153137208, -93.7744873626709, -276.5351152771,
    82.2280757629395, -413.176540896606, -90.4916240661622, -105.857697381592,
    154.842738793945, 849.984451940918, 101.827337915039, 1874.18971590271,
    1368.88667627563, -57.3683639556885, 1275.62188513184, -89.2440485870361,
    -417.936036978149, 392.281618716431, 55.5831157928466, -41.3011304901123,
    206.429591436768, 28.1189699615479, -355.684031378174, -11.9870409210205,
    360.547564938355, -1739.89044494934, 192.427370782471, -2.51351051330565,
    -349.780997296143, -10.8496227722167, 237.020107214355, 805.973101419068,
    751.231902365112, -0.764981460571292, -311.423515466309, 108.36333951416,
    216.091088717651, 8.15076817932131, 1464.24882801819, -8.44452106018079,
    -228.265222247314, 171.337487612915, -46.7040852630615, -252.735012075806,
    1739.34752667847, 69.7715547912598, -116.013154119873, -199.826271469116,
    163.103664303589, 137.761358477783, -121.910942614746, -139.66812940979,
    25.0441816223144, 161.536108007812, 9.45779364624022, -230.24455713501,
    -27.0637326782225, -60.4999794891359, 21.8207683135986, -121.908319821167,
    215.133769061279, 356.064336447144, -56.2746590332031, 36.9691497619629,
    285.947447167969, 206.478550250244, -1590.6648557373, -297.025252981567,
    1342.76365222778, -412.473632217407, -446.423072305298, -243.027178775024,
    -450.760298620605, 329.826783746338, -143.601445513916, -88.5385171142579,
    -1668.706952948, -11.1014109558106, -7.66642563171389, -1562.65779163513,
    -227.605152529907, -1139.62741526184, 245.455011364746, 359.272887258911,
    120.45878923645, 201.233837356567, 466.392148352051, 346.973733901978,
    237.866395275879, -1356.37070531616, 36.7593262756348, -12.1216776580811,
    -765.456186209107, -1432.98775134888, 97.4516439605713, -21.200040499878,
    -1066.9917698822, -553.989082571411, -82.8426837249756, 1077.79330810547,
    -375.020139907837, -366.178702752686, -528.385371652222, 275.786744842529,
    68.3115330322268, -421.008202523804, -938.025512539673, -403.418874517822,
    -1273.43447528687, -137.444874719238, -1160.50572641602, 236.189555914307,
    235.681608224487, 400.36069720459, 482.495226663208, 363.804200299072,
    427.917515075684, -1472.14343669128, 105.337509988403, 221.325310437012,
    -1599.68901417847, -391.493032113648, -271.700432446289, -1407.48982643738,
    0.023900994106419, 0.0231330567647606, 30.2897002641407, 59.9912724843918,
    16.1416492585762, 30.0108257396054, 7.73999646008463, 8.6370957784742,
    2.52437179291348, 3.25387945289847, 3.22067975233017, 1.48062419886981,
    5.06134927258143, 3.17062128567741, 3.36351305458512, 4.80438388980674,
    3.17258160349083, 47.2580790486122, 0.0289358944731071, 0.0268201565710897,
    0.0225909765801275, 0.0260099506467678, 0.0907407608883197,
    0.0761823890648398, 8.98772516886864, 28.5410912655296, 6.79983236422023,
    6.49731673520026, 9.81290346155874, 2.11201634584801, 27.287196596754,
    12.6062100672994, 0.922579033608156, 0.766831146015318, 0.246310808643451,
    0.157133320818368, 4.50291965180116, 30.5014231378917, 70.4029586552732,
    2.02927675028875, 0.585308131393753, 0.175132418856374, 0.274965181689283,
    0.123785091963862, 0.166903875432994, 0.213032216677792, 0.240248928089632,
    0.1187406225148, 0.129183503389083, 0.715304161233158, 1.50504390388076,
    0.248707421092679, 0.195528418775028, 2.31111769071423, 1.60878616349495,
    0.170860037069368, 0.306429009677705, 0.372973925270621, 3.89355640648065,
    0.849794129491282, 3.2635762050366, 3.28394909256719, 1.42194806139652,
    3.47618472716798, 2.01514804679201, 0.727915104532024, 0.964802010654108,
    0.806483183851841, 0.623699982478079, 2.25441017209958, 0.610253465224155,
    1.57185804485401, 0.803193766967432, 2.56177702490827, 6.27705719612543,
    7.27302055605245, 12.3025469039091, 3.19083429783221, 2.55035834402969,
    5.66383214748435, 2.26441302819176, 48.7185560480955, 38.064614949065,
    31.1043982524214, 7.89212651686726, 3.32664849130564, 6.42201175187935,
    1.46025487631199, 10.6826245135661, 43.1109563076373, 58.1062991926117,
    42.2165131685901, 38.5909605614048, 36.5339163848374, 39.0250044607596,
    37.389874638775, 34.8185396852934, 37.2046832077029, 38.5269897435591,
    52.0285431329025, 37.0800207785229, 38.2543736124358, 40.6883215556404,
    44.9868977425105, 41.6355896555399, 23.424977946178, 37.5206596857773,
    46.2511382671431, 24.325899069968, 22.8365819537987, 35.9081658546592,
    40.5635359156785, 49.7031914392606, 41.8707686221091, 37.3920591291155,
    39.246866094806, 36.9839188943937, 47.7222387407682, 43.5548016056201,
    47.6661850701958, 48.2394849129995, 49.3781948348789, 41.3575416927592,
    44.1770071819466, 39.2826066302928, 36.0873421343397, 43.5256706344402,
    39.7814615371494, 38.6567023952772, 36.6978695946594, 44.7160166930849,
    42.0350978026939, 60.7290535961573, 58.6884511493973, 47.0441544038059,
    52.0429494195169, 40.2325920407092, 39.8428189850556, 44.7757742626577,
    45.4614369771841, 44.1189123848991, 42.1823325832929, 42.966645439567,
    40.0754814251924, 46.0118737075095, 37.226456950665, 37.6400898427513,
    40.7161692803813, 41.3009481887216, 38.6207180738967, 44.0325801970388,
    43.6332008349993, 40.2469424522467, 44.528264981689, 43.4040539063303,
    44.884410053672, 43.9303093008641, 41.0205892113336, 40.7678160757023,
    50.5213004528575, 40.0620620426386, 39.4044365394214, 38.6079443863076,
    40.012189361057, 39.3708833677596, 53.9186504243407, 41.4836564404552,
    44.5602756392325, 42.1733287715947, 42.0816958035646, 43.0055716825984,
    37.5721112836435, 38.7317086791101, 37.9477339208664, 38.6649543737248,
    38.597272263127, 38.7985162397143, 39.6648748098574, 41.0540518879328,
    41.9856036638906, 42.8627844572998, 41.5986959111667, 37.3472008270586,
    46.2602777467611, 40.8229131039719, 42.8330205862047, 40.6315551426751,
    53.0944117952853, 66.9743357469595, 75.4056194186102, 52.9049016280758,
    52.6597445665472, 56.1866414286977, 50.0491625151176, 62.6259867558734,
    60.6372516500431, 53.1266888601478, 62.8801450923029, 48.3432761028053,
    90.4771965230962, 83.1949024312402, 48.3892788378523, 80.944087597953,
    62.935997879368, 61.2508468330903, 56.2933400092884, 55.1306388788584,
    54.5730939815481, 54.7043814234442, 54.2382262849893, 66.9631234531142,
    55.7866155516914, 48.0236416231127, 59.3138970680715, 65.0847904382007,
    53.0623488240625, 54.8521633795295, 56.2900965515109, 54.6566038499067,
    56.3193210625947, 83.2665554001824, 52.8711005659689, 50.3456559590632,
    51.2394650577083, 65.1449055255141, 56.4180151979781, 53.593422491879,
    56.8324611867939, 66.3305213683673, 78.3665174605801, 44.3309688863162,
    64.3505032669697, 68.0314021702132, 72.1742005650908, 64.6382973270138,
    57.7878574147254, 57.4874773613056, 52.6560801784698, 56.5380804993766,
    53.0611439866183, 51.7593311978994, 61.2032355348196, 107.554063069431,
    82.6540449398095, 126.691306601288, 0.0393419036865234, 0.0393419036865234,
    179.470770501709, 144.354187271118, 66.2342805175781, 103.404511120605,
    12.661973135376, 13.8518471557617, 6.85161109313965, 12.171510736084,
    5.84795541687012, 4.7052916809082, 16.7211833312988, 11.4021579528809,
    10.4920485809326, 10.9116955535889, 6.50190528259277, 148.204448245239,
    0.0629470458984375, 0.0603242523193359, 0.0603242523193359,
    0.0603242523193359, 0.340088900756836, 0.270147738647461, 46.9821013824463,
    68.49600284729, 15.0390983825684, 16.7194348022461, 19.4952246734619,
    3.58186176452637, 73.0684063201904, 32.6083183044434, 4.09592930603027,
    4.4229042388916, 0.597122671508789, 0.457240347290039, 8.59839161682129,
    119.665831311035, 41.6158657196045, 7.38578671875, 2.22937454223633,
    0.457240347290039, 0.92322333984375, 0.340088900756836, 0.456366082763672,
    0.549912387084961, 0.736130731201172, 0.340963165283203, 0.386424920654297,
    0.433635205078125, 2.13495397338867, 0.5035763671875, 0.502702102661133,
    7.47845875854492, 2.92878616333008, 0.479971224975586, 0.948577011108398,
    1.2283416595459, 14.2478889862061, 4.35383734130859, 5.59092164611816,
    5.07772836914063, 2.51176198425293, 12.4276702423096, 5.12406438903809,
    2.1847870513916, 3.6316948425293, 2.65077004394531, 2.37100539550781,
    11.0970396331787, 2.4645516998291, 7.31759408569336, 2.27745909118652,
    9.6746112487793, 17.2807126281738, 21.900326385498, 46.6577492431641,
    15.7647379394531, 6.61818246459961, 11.0270984710693, 5.14766953125,
    163.207701782227, 99.717737612915, 83.2439711425781, 12.4040651000977,
    8.13415715332031, 14.7138719787598, 3.88785434875488, 26.8713944824219,
    72.0734932891846, 169.186796878052, 84.0902592041016, 66.4274929779053,
    61.643517489624, 54.3871219207764, 66.3566775512695, 69.2032828491211,
    60.6398618133545, 68.5965432678223, 163.586258322144, 69.7636864105225,
    63.6263494354248, 76.0409057098389, 108.333614520264, 118.436615386963,
    50.8874410217285, 93.1896043945312, 147.977139468384, 48.9509450958252,
    53.3135250823975, 59.007609942627, 80.2399982299805, 174.179721588135,
    70.8836192687988, 60.8269544219971, 63.4174002136231, 60.3837023071289,
    85.2574023468018, 62.6707783081055, 75.3642249664307, 98.579445199585,
    120.793632550049, 64.46739190979, 69.6001989440918, 65.8207533966065,
    65.1895344085693, 85.0939148803711, 52.0773150421143, 61.2238705169678,
    69.2504931335449, 107.796816101074, 61.4336940032959, 177.166209210205,
    163.236552511597, 103.129992059326, 141.839802493286, 65.843484274292,
    52.0064996154785, 84.3472929748535, 77.0200819793701, 85.2801332244873,
    71.560300012207, 67.8971316467285, 65.283954977417, 84.3464187103271,
    69.8799635925293, 36.9577843231201, 76.6240401489258, 70.2532745452881,
    76.1099726074219, 67.7336441802979, 72.2841910400391, 92.0696715362549,
    109.429942236328, 67.7563750579834, 64.6999462738037, 64.3738456054687,
    79.634132913208, 73.426854776001, 154.907434368896, 60.3373662872315,
    67.9906779510498, 84.4399650146484, 72.8900563568115, 63.4637362335205,
    133.510684350586, 62.6471731658936, 61.9241564025879, 54.9466512176514,
    63.6272236999512, 72.9599975189209, 75.1998632354736, 66.8008039306641,
    71.3968125457764, 61.4572991455078, 59.0539459625244, 69.8808378570557,
    59.8469038879395, 55.7632142852783, 73.8010399932861, 60.8269544219971,
    61.0367779083252, 67.4075435119629, 71.4431485656738, 68.4103249237061,
    72.1207035736084, 73.5903422424317, 58.7033658874512, 106.543120770264,
    221.294711178589, 94.129438760376, 110.229020013428, 99.5883464630127,
    81.2252943511963, 113.985734683228, 96.8352874694824, 98.3757415649414,
    115.501709371948, 97.4656321929932, 285.01810397644, 146.511872122192,
    106.752069992065, 113.822247216797, 129.78194614563, 147.935174771118,
    103.811918389893, 112.77225552063, 102.108851092529, 113.169171615601,
    124.322164178467, 88.0157069274902, 116.29554156189, 93.638976361084,
    85.5887486022949, 85.5424125823975, 109.995591384888, 101.922632748413,
    93.381942590332, 81.9255802368164, 94.2457159423828, 212.894777609253,
    84.6323032104492, 88.9249420349121, 81.3188406555176, 135.358879559326,
    118.232037487793, 131.648500909424, 86.6859505828857, 229.695519012451,
    163.171856936646, 95.5754722869873, 96.0886655639648, 123.085079873657,
    150.968872677612, 140.655174060059, 123.132290158081, 124.088735549927,
    121.848869833374, 117.555356744385, 118.559012420654, 126.818189401245,
    84.4688157440186, 214.17819793396, 153.418561880493, 244.908596035767,
    -0.0542044006347656, -0.0542044006347656, -3.85725509033203,
    -214.159838378906, -45.6007634307861, -90.143666784668, -37.317107043457,
    -41.5406789703369, -9.64401199035645, -6.5412471862793, -11.1136506591797,
    -3.15696920471191, -18.6742902832031, -15.2672814239502, -7.54402859802246,
    -28.940778616333, -17.8804580932617, -217.099989981079, -0.100540420532227,
    -0.0542044006347656, -0.0340963165283203, -0.0340963165283203,
    -0.31386096496582, -0.312986700439453, -7.07979413452148, -143.066395623779,
    -17.4398287719727, -15.3162402374268, -38.8601839324951, -6.05340758056641,
    -98.5234922698975, -39.0464022766113, -1.01327258605957, -0.570020471191406,
    -0.687171917724609, -0.196709518432617, -22.4563586242676, -33.8663849578857,
    -326.162741116333, -3.51017207336426, -1.33937325439453, -0.406533004760742,
    -0.663566775512695, -0.220314660644531, -0.267524945068359,
    -0.407407269287109, -0.383802127075195, -0.196709518432617,
    -0.266650680541992, -2.48378551940918, -7.03345811462402, -0.45374328918457,
    -0.569146206665039, -5.32951655273438, -4.37307116088867, -0.336591842651367,
    -1.0578600769043, -0.707280001831055, -7.24065880737305, -2.03791061096191,
    -12.4914915527344, -16.5743068908691, -6.33142369995117, -6.42497000427246,
    -7.35781025390625, -1.36122986755371, -1.31489384765625, -2.31767525939941,
    -1.24407842102051, -5.95723848266602, -1.26768356323242, -4.55754097595215,
    -3.5783647064209, -8.7819871673584, -18.4181307769775, -26.3039968048096,
    -30.1079217590332, -1.87442314453125, -8.75838202514649, -15.4080380126953,
    -7.80106236877441, -223.12104977417, -102.83798770752, -86.7139270477295,
    -37.014611517334, -12.74765105896, -20.5181141693115, -3.53115442199707,
    -42.8477044372559, -75.4385374511719, -82.8348153442383, -76.5112600250244,
    -77.9118317962647, -83.4415549255371, -81.2480252288818, -81.8075545257568,
    -79.5440836669922, -84.6777649658203, -84.9347987365723, -49.514845715332,
    -60.2945273254395, -82.8811513641358, -85.7522360687256, -62.4417210021973,
    -57.1917625213623, -53.9246359863281, -43.540996206665, -47.6019549316406,
    -53.1080729187012, -39.9713741455078, -85.9620595550537, -74.107906842041,
    -43.4247190246582, -74.7609824432373, -68.624519732666, -88.0847738250732,
    -69.0913769897461, -83.9547482025147, -76.441318862915, -94.8043709747314,
    -82.9519667907715, -77.6311928833008, -79.6848402557373, -93.5681609344482,
    -76.9309069976807, -69.0214358276367, -98.46753934021, -84.841252432251,
    -81.2716303710938, -72.2876880981445, -61.2046366973877, -81.1081429046631,
    -76.7210835113525, -64.1919985839844, -81.3415715332031, -91.7479421905518,
    -87.3853622039795, -105.771145193481, -93.2647911437988, -88.1083789672852,
    -100.125144882202, -101.384960064697, -85.6350846221924, -85.565143460083,
    -105.702078295898, -92.2847406097412, -124.647390582275, -82.3216220672608,
    -88.948547177124, -85.6114794799805, -91.117597467041, -85.8213029663086,
    -72.3812344024658, -75.7419072418213, -81.2016892089844, -101.688329855347,
    -86.8013535003662, -86.4044374053955, -84.2808488708496, -63.2346789276123,
    -80.7112268096924, -89.3218581298828, -73.057915145874, -75.2741757202148,
    -78.0517141204834, -84.1409665466309, -115.151129296875, -109.107338626099,
    -106.774800869751, -89.205580947876, -97.9779512054443, -73.2214026123047,
    -75.5784197753906, -72.1486800384521, -74.3177303283691, -86.7777483581543,
    -72.8716968017578, -74.598369241333, -101.618388693237, -87.2918158996582,
    -94.6644886505127, -101.571178408813, -75.6710918151856, -98.958001739502,
    -81.2244200866699, -99.191430368042, -70.1877047058106, -178.547547161865,
    -157.845837442017, -123.172506326294, -124.339649468994, -115.355707196045,
    -123.498606994629, -96.8527727600098, -103.035571490479, -146.178777337647,
    -130.755876828003, -163.959569274902, -84.4163598724365, -182.275411102295,
    -252.205207772827, -110.829639743042, -226.16523885498, -112.392824716187,
    -105.135554882813, -111.692538830566, -100.959193240356, -98.9990921722412,
    -86.1657631896973, -96.5257978271484, -250.315047866821, -96.712890435791,
    -102.219882687378, -174.109780426025, -208.245438858032, -98.6030503417969,
    -111.459110202026, -169.559233566284, -141.162247485352, -146.389475088501,
    -124.899178765869, -97.2024785705566, -110.129353857422, -132.062028030396,
    -105.019277700806, -124.269708306885, -114.843388183594, -171.355847167969,
    -166.152224707031, -189.321983184814, -75.6894513702393, -178.96282281189,
    -98.6965966461182, -104.2962609375, -130.569658483887, -93.0261169281006,
    -83.0621241210938, -86.4927381225586, -185.355445028687, -94.8227305297852,
    -73.986384072876, -154.718593231201, -273.438470324707, -143.449323486328,
    -278.805580252075, -0.0533301361083984, -0.171355847167969, 312.911513690186,
    -475.162770080566, 41.61936277771, 68.9987049499512, -101.801109979248,
    38.1883116439819, -13.3753729888916, 33.7571018920899, 26.7651713424683,
    -16.3697289916992, -61.1985168457031, 3.96959808197022, 26.9723720352173,
    -14.4839404083252, -6.63916481323242, -91.5267532653809, 0.195398121643066,
    0.233865760803223, 0.264902151489258, 0.383802127075195, 0.360196984863281,
    -0.155619085693359, 75.3860815795898, -132.178742344666, 118.602725646973,
    48.3293430175781, -145.909940995789, -1.11337587432861, 135.890432391357,
    -111.146123501587, 30.1433294723511, 12.4075621582031, 0.500516441345215,
    4.42071857757568, -20.4674068267822, 775.5473845047, -1070.30654383392,
    33.5477155380249, 3.87954883575439, -1.30877399597168, 0.404784475708008,
    1.79136801452637, 1.80404485015869, -0.354951397705078, 3.56262794494629,
    1.91201651916504, 2.8216887588501, -10.3294353790283, -14.5630613479614,
    1.22615599822998, 5.17564599609375, 10.971582673645, 9.14917826843262,
    2.7390707611084, -1.29697142486572, 2.4807255935669, 37.9841708770752,
    -3.94118448486328, -33.1210744491577, -31.9163379318237, -0.128516885375977,
    3.74010364379883, -8.74045960235596, 6.37513692626953, 1.37784089355469,
    0.715585514831543, -2.46892302246094, 7.69658775787354, 6.1176660232544,
    1.86699189605713, 1.50985483703613, 14.8415145996094, -39.9058043060303,
    2.70934576721191, -52.1651786270142, 38.8142850448608, -12.8320175857544,
    3.81660178985595, -26.2952541595459, -148.782337097168, 108.947785350037,
    96.6971536743164, 35.5065052093506, 63.0593888900757, 74.8138754470825,
    24.2127560577393, -63.2748950958252, 231.262638175964, 1278.59045033112,
    -180.150948303223, 123.414240467835, 166.940374177551, 114.024639454651,
    -333.034460293579, 305.908654833984, -289.84579269104, 339.875580212402,
    1193.95508719482, 757.389347424316, -214.978587107849, 293.314874331665,
    735.455361854553, 601.601091545105, -31.8848644088746, 1180.55610906372,
    1035.46841385498, -2.7093457672119, -0.817874464416512, -799.033626741028,
    -334.194609320068, 2930.59414237061, 276.937714091492, -239.517880966187,
    -3.87124332275388, 155.800495582581, -378.253607258606, 389.76198835144,
    -439.531682176209, 589.755244345093, 706.326616365051, -13.3989781311035,
    -81.4635314346314, -282.295644241333, 107.31291068573, -414.329258674622,
    -94.817484942627, -112.680457745361, 120.650690299988, 821.100063386536,
    112.808100366211, 1858.18542948303, 1322.91085336304, -71.5874022125244,
    1278.72770986176, -53.4490360839843, -443.68137875061, 395.200350837708,
    74.4737865463257, -56.9670765380859, 235.294309039307, 30.9677609207153,
    -366.473329898071, -27.0230793777466, 345.640480499268, -1743.25986043396,
    197.154519076538, -23.6982513839722, -368.736363624573, -2.57558329467773,
    224.749367454529, 773.904204327393, 744.794255525207, 9.84247003784179,
    -323.857742692566, 136.272048857117, 216.769955122376, -3.0935850265503,
    1439.3895533432, 5.38240955657957, -253.730799371338, 201.194495452881,
    -49.4785637374878, -282.295644241333, 1749.61620067291, 93.958082913208,
    -157.602354771423, -228.266970776367, 154.496092909241, 153.968474267578,
    -116.671038175964, -172.05219886322, 34.5898388534546, 204.494406907654,
    9.22611354675289, -256.211962097168, -9.77165461120603, -52.7164024108887,
    32.9296105178833, -153.382279902649, 210.807908184814, 386.620755908203,
    -47.2190270690918, 3.64000035552976, 323.32968691864, 213.609051727295,
    -1612.75839458313, -318.915088192749, 1327.66379246063, -363.31373789978,
    -431.446483836365, -231.620649499512, -447.088387609863, 354.463558099365,
    -141.726148104858, -161.434256190491, -1695.95165838318, -55.5306599212647,
    -3.87823743896489, -1546.98091728058, -185.293372247314, -1104.03523212891,
    163.657510881042, 301.450780014038, 64.6199510696411, 221.273728829956,
    474.415711042786, 383.521925294495, 254.593261326599, -1353.8768657547,
    46.2529647674561, 4.59338582153316, -720.181086314392, -1413.08774219971,
    118.331703643799, -32.1344669311524, -1049.79848370667, -539.315863893128,
    -58.2644850952148, 1117.60818890076, -319.630673707581, -333.352255448914,
    -481.955368318176, 346.739431008911, 40.7166846542358, -414.02282895813,
    -987.401350195313, -360.535762367249, -1287.27451987152, -118.256079762268,
    -1160.50878634186, 249.591593971252, 286.521401829529, 381.247963261414,
    430.121535946655, 338.391516178894, 460.674458349609, -1415.89063427124,
    118.972976673889, 181.492944351196, -1514.66635185699, -301.558314550781,
    -247.541880789185, -1346.35993935699, 0.0, 0.0, 1.0, 2.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 3.0,
    3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 6.0, 5.0, 5.0, 3.0, 2.0, 2.0, 2.0, 3.0,
    6.0, 5.0, 6.0, 5.0, 5.0, 4.0, 3.0, 0.0, 7.0, 3.0, 0.0, 1.0, 2.0, 2.0, 7.0,
    4.0, 1.0, 2.0, 3.0, 3.0, 4.0, 5.0, 4.0, 2.0, 4.0, 5.0, 3.0, 2.0, 3.0, 2.0,
    2.0, 3.0, 5.0, 5.0, 5.0, 6.0, 4.0, 4.0, 5.0, 2.0, 4.0, 4.0, 4.0, 8.0, 5.0,
    3.0, 3.0, 4.0, 0.0, 4.0, 6.0, 2.0, 6.0, 6.0, 5.0, 7.0, 7.0, 4.0, 5.0, 5.0,
    5.0, 3.0, 4.0, 1.0, 5.0, 3.0, 4.0, 7.0, 2.0, 1.0, 3.0, 4.0, 4.0, 2.0, 4.0,
    4.0, 5.0, 3.0, 4.0, 6.0, 3.0, 3.0, 4.0, 3.0, 5.0, 5.0, 3.0, 5.0, 5.0, 1.0,
    4.0, 5.0, 6.0, 7.0, 5.0, 5.0, 8.0, 7.0, 6.0, 2.0, 6.0, 3.0, 6.0, 5.0, 6.0,
    7.0, 8.0, 9.0, 7.0, 9.0, 7.0, 7.0, 2.0, 8.0, 7.0, 5.0, 5.0, 8.0, 7.0, 4.0,
    5.0, 7.0, 6.0, 8.0, 4.0, 6.0, 9.0, 8.0, 6.0, 4.0, 7.0, 5.0, 5.0, 3.0, 9.0,
    8.0, 8.0, 10.0, 7.0, 9.0, 3.0, 7.0, 7.0, 4.0, 8.0, 8.0, 8.0, 0.0, 0.0, 0.0,
    4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 9.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 8.0, 19.0, 0.0, 0.0, 0.0, 0.0, 0.0, 36.0, 36.0, 36.0,
    37.0, 20.0, 16.0, 20.0, 40.0, 25.0, 43.0, 20.0, 40.0, 24.0, 42.0, 25.0, 40.0,
    0.0, 30.0, 8.0, 0.0, 0.0, 22.0, 17.0, 42.0, 38.0, 0.0, 16.0, 35.0, 20.0,
    42.0, 23.0, 42.0, 3.0, 36.0, 38.0, 19.0, 2.0, 22.0, 2.0, 20.0, 5.0, 21.0,
    38.0, 30.0, 43.0, 24.0, 33.0, 39.0, 2.0, 42.0, 39.0, 36.0, 37.0, 37.0, 20.0,
    36.0, 40.0, 0.0, 39.0, 40.0, 19.0, 37.0, 40.0, 25.0, 41.0, 39.0, 18.0, 36.0,
    37.0, 23.0, 18.0, 37.0, 0.0, 40.0, 38.0, 19.0, 38.0, 4.0, 0.0, 36.0, 38.0,
    38.0, 18.0, 20.0, 36.0, 37.0, 36.0, 18.0, 34.0, 41.0, 20.0, 21.0, 36.0, 22.0,
    37.0, 17.0, 41.0, 39.0, 0.0, 29.0, 25.0, 42.0, 41.0, 41.0, 34.0, 41.0, 39.0,
    38.0, 11.0, 37.0, 38.0, 46.0, 38.0, 26.0, 29.0, 38.0, 39.0, 39.0, 42.0, 40.0,
    37.0, 13.0, 39.0, 41.0, 29.0, 28.0, 41.0, 38.0, 30.0, 25.0, 41.0, 40.0, 46.0,
    40.0, 39.0, 40.0, 39.0, 36.0, 22.0, 40.0, 46.0, 41.0, 12.0, 40.0, 40.0, 39.0,
    45.0, 38.0, 39.0, 36.0, 37.0, 47.0, 36.0, 43.0, 36.0, 40.0, 0.0, 0.0, 0.0,
    2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 2.66666666666667, 6.33333333333333, 0.0, 0.0, 0.0, 0.0,
    0.0, 7.2, 6.0, 7.2, 7.4, 6.66666666666667, 8.0, 10.0, 20.0, 8.33333333333333,
    7.16666666666667, 4.0, 6.66666666666667, 4.8, 8.4, 6.25, 13.3333333333333,
    0.0, 4.28571428571429, 2.66666666666667, 0.0, 0.0, 11.0, 8.5, 6.0, 9.5, 0.0,
    8.0, 11.6666666666667, 6.66666666666667, 10.5, 4.6, 10.5, 1.5, 9.0, 7.6,
    6.33333333333333, 1.0, 7.33333333333333, 1.0, 10.0, 1.66666666666667, 4.2,
    7.6, 6.0, 7.16666666666667, 6.0, 8.25, 7.8, 1.0, 10.5, 9.75, 9.0, 4.625, 7.4,
    6.66666666666667, 12.0, 10.0, 0.0, 9.75, 6.66666666666667, 9.5,
    6.16666666666667, 6.66666666666667, 5.0, 5.85714285714286, 5.57142857142857,
    4.5, 7.2, 7.4, 4.6, 6.0, 9.25, 0.0, 8.0, 12.6666666666667, 4.75,
    5.42857142857143, 2.0, 0.0, 12.0, 9.5, 9.5, 9.0, 5.0, 9.0, 7.4, 12.0, 4.5,
    5.66666666666667, 13.6666666666667, 6.66666666666667, 5.25, 12.0, 4.4, 7.4,
    5.66666666666667, 8.2, 7.8, 0.0, 7.25, 5.0, 7.0, 5.85714285714286, 8.2, 6.8,
    5.125, 5.57142857142857, 6.33333333333333, 5.5, 6.16666666666667,
    12.6666666666667, 7.66666666666667, 7.6, 4.33333333333333, 4.14285714285714,
    4.75, 4.33333333333333, 5.57142857142857, 4.66666666666667, 5.71428571428571,
    5.28571428571429, 6.5, 4.875, 5.85714285714286, 5.8, 5.6, 5.125,
    5.42857142857143, 7.5, 5.0, 5.85714285714286, 6.66666666666667, 5.75, 10.0,
    6.5, 4.44444444444445, 4.875, 6.0, 5.5, 5.71428571428571, 9.2, 8.2, 4.0,
    4.44444444444445, 5.0, 4.875, 4.5, 5.42857142857143, 4.33333333333333, 12.0,
    5.28571428571429, 6.71428571428571, 9.0, 5.375, 4.5, 5.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0,
    4.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0, 0.0, 0.0, 0.0, 3.0,
    4.0, 0.0, 1.0, 5.0, 0.0, 0.0, 0.0, 0.0, 14.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 6.0, 6.0,
    1.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 7.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 8.0, 1.0, 1.0, 1.0, 0.0,
    4.0, 1.0, 2.0, 1.0, 1.0, 2.0, 2.0, 1.0, 1.0, 4.0, 4.0, 3.0, 3.0, 3.0, 3.0,
    3.0, 0.0, 2.0, 1.0, 0.0, 0.0, 2.0, 4.0, 1.0, 0.0, 1.0, 6.0, 0.0, 0.0, 0.0,
    4.0, 2.0, 3.0, 0.0, 0.0, 4.0, 1.0, 1.0, 8.0, 8.0, 5.0, 4.0, 6.0, 3.0, 1.0,
    3.0, 3.0, 0.0, 7.0, 9.0, 4.0, 50.0, 50.0, 48.0, 49.0, 50.0, 49.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 49.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 49.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 48.0, 49.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 46.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    44.0, 50.0, 50.0, 50.0, 47.0, 46.0, 50.0, 49.0, 45.0, 50.0, 50.0, 50.0, 50.0,
    36.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 49.0, 46.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 47.0, 50.0, 43.0, 44.0, 49.0, 46.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 49.0, 48.0, 50.0, 50.0, 50.0, 50.0, 50.0, 44.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 43.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 48.0, 41.0,
    49.0, 49.0, 49.0, 50.0, 46.0, 49.0, 48.0, 49.0, 49.0, 46.0, 48.0, 49.0, 49.0,
    46.0, 46.0, 47.0, 47.0, 47.0, 47.0, 47.0, 50.0, 48.0, 49.0, 50.0, 50.0, 48.0,
    46.0, 49.0, 50.0, 49.0, 41.0, 50.0, 50.0, 50.0, 46.0, 48.0, 47.0, 50.0, 49.0,
    46.0, 49.0, 49.0, 42.0, 42.0, 45.0, 46.0, 44.0, 47.0, 49.0, 47.0, 47.0, 50.0,
    40.0, 41.0, 42.0, 0.0236051422119141, 0.000874264526367186, 29.585985836792,
    110.7387162323, 3.94205874938965, 11.4563623535156, 7.72237856140137,
    2.09998339233398, 1.09720198059082, 1.37696662902832, 2.40422744750977,
    0.466857257080078, 4.38705939331055, 1.14353800048828, 2.75305899353027,
    0.838419680786132, 9.61166420288086, 6.06652154846191, 0.0218566131591797,
    0.0227308776855469, 0.0454617553710937, 0.0454617553710937,
    0.0227308776855469, 0.0227308776855469, 6.27634503479004, 5.11007615661621,
    6.13646271057129, 0.489588134765622, 0.723891027832028, 0.0944205688476556,
    7.23191616210938, 5.01740411682129, 0.208949221801758, 0.980050534057617,
    0.13988232421875, 0.069941162109375, 15.5400519561768, 24.7102125732422,
    21.8155227264404, 1.84294962158203, 0.0454617553710937, 0.0708154266357423,
    0.0227308776855469, 0.0690668975830079, 0.0926720397949219,
    0.0236051422119141, 0.0944205688476562, 0.047210284423828, 0.162613201904297,
    0.326974932861328, 1.28342032470703, 0.0463360198974609, 0.0463360198974609,
    1.05086596069336, 1.46963866882324, 0.116277182006836, 0.0935463043212892,
    1.42242838439941, 2.70584870910644, 0.442377850341797, 4.57240347290039,
    0.79295792541504, 5.78675690002442, 1.39969750671387, 0.209823486328125,
    0.13988232421875, 0.465982992553711, 0.208949221801758, 0.186218344116211,
    3.8021764251709, 0.2797646484375, 1.18987402038574, 0.0944205688476563,
    0.63121898803711, 3.40613459472656, 0.0699411621093757, 1.98283194580078,
    0.209823486328125, 1.04999169616699, 2.49689948730469, 0.676680743408203,
    22.1669770660401, 10.1030008666992, 95.3158157226563, 0.630344723510741,
    0.325226403808594, 0.466857257080078, 2.07637825012207, 12.7153032714844,
    4.03647931823731, 5.0864710144043, 3.3597985748291, 1.51684895324707,
    0.232554364013673, 3.14997508850097, 1.16714314270019, 8.07208437194824,
    24.4313221893311, 9.72968991394043, 8.16475641174316, 7.27825218200683,
    6.55610968322755, 6.9302949005127, 10.5698581237793, 1.35511001586914,
    4.52606745300292, 5.53059739379883, 8.56516956481934, 3.71037864990233,
    6.37163986816406, 7.95668145446778, 9.14655547485352, 3.22079051513671,
    9.89317738037109, 1.47051293334962, 9.84684136047363, 1.16539461364746,
    1.75027758178711, 8.19011008300782, 2.07637825012207, 20.3450097930908,
    16.5437076324463, 8.0974380432129, 15.5155725494385, 5.18089158325195,
    16.6827156921387, 10.96677421875, 3.45334487915039, 11.9459504882813,
    1.63400039978027, 12.2974048278809, 2.16992455444336, 8.51533648681641,
    0.187092608642573, 8.4462695892334, 4.94746295471192, 4.92473207702638,
    5.1581607055664, 1.56318497314453, 3.26712653503418, 1.51684895324706,
    3.70950438537599, 7.04657208251953, 3.71125291442871, 11.6670601043701,
    11.526303515625, 4.6440931640625, 4.59775714416504, 3.45247061462403,
    5.71681573791505, 3.21991625061035, 6.83762286071777, 7.69877341918946,
    3.68764777221681, 12.9968164489746, 1.91289078369141, 3.08090819091797,
    3.12724421081543, 3.05555451965332, 18.5746241271973, 1.93824445495605,
    6.44070676574707, 1.98283194580079, 5.18001731872558, 4.29351308898926,
    4.12915135803223, 2.77753840026855, 16.519228225708, 10.3355552307129,
    16.8007414031982, 35.0230369262695, 0.956445391845698, 5.50699225158691,
    0.467731521606446, 2.14544514770508, 1.42330264892578, 0.840168209838872,
    20.2296068756104, 3.91932787170411, 8.07295863647461, 8.56429530029296,
    7.41900877075196, 8.14289979858398, 4.47885716857911, 3.82665583190919,
    1.84294962158204, 0.34970581054688, 3.64043748779297, 29.516918939209,
    38.0794657104492, 42.2803067596435, 33.4834570953369, 7.07017722473145,
    22.2369182281494, 78.2134530578613, 17.335791293335, 25.3851447875977,
    0.885629965209944, 51.2170387481689, 43.2830881713868, 36.3291881286621,
    6.23088327941895, 45.2895252593994, 25.6421785583496, 10.3836397796631,
    14.3965139556885, 18.9933968353271, 28.1163471679687, 13.6035560302735,
    8.63511072692873, 43.4229704956055, 17.3139346801758, 16.77626199646,
    20.719195010376, 35.793263973999, 16.4973716125488, 38.1266759948731,
    0.186218344116213, 25.2688676055908, 57.4925095184326, 1.04999169616698,
    45.2423149749756, 29.2371542907715, 48.7909546875, 3.73310952758789,
    41.8143237670898, 7.25726983337404, 1.14353800048826, 13.1603039154053,
    23.2624305175781, 13.2765810974121, 15.7498754425049, 56.1155428894043,
    1.68033641967773, 4.26990794677738, 29.773952709961, 14.0232030029297,
    42.0460038665771, 44.2168026855469, 0.162613201904293, 8.39993356933596,
    46.7836433349609, 95.8534884063721, 41.4392642852783, 78.9828058410645,
    0.00124237590589021, 0.000124894932338169, 29.585985836792, -27.684679058075,
    0.303235288414589, -3.81878745117188, 3.86118928070069, 2.09998339233398,
    0.274300495147705, 1.37696662902832, -0.400704574584961, -0.233428628540039,
    2.19352969665528, 1.14353800048828, 2.75305899353027, -0.838419680786132,
    -1.60194403381348, 6.06652154846191, -0.00546415328979492,
    0.000909235107421875, 0.00239272396689967, 0.00168376871744792,
    0.0227308776855469, 0.0227308776855469, -6.27634503479004, 2.5550380783081,
    -3.06823135528564, -0.489588134765622, -0.0329041376287286,
    -0.0472102844238278, 1.44638323242188, -0.135605516670846,
    0.0189953838001598, -0.980050534057617, 0.13988232421875, 0.0116568603515625,
    0.914120703304515, -12.3551062866211, 7.27184090881348, -1.84294962158203,
    0.0227308776855468, 0.0354077133178711, 0.00133711045209099,
    -0.00181754993639494, -0.0926720397949219, 0.0236051422119141,
    0.00219582718250363, -0.047210284423828, -0.00416956927959736,
    -0.163487466430664, -0.641710162353516, 0.00463360198974609,
    -0.0115840049743652, -0.525432980346679, 0.244939778137207,
    -0.00322992172241211, 0.0935463043212892, -0.711214192199707,
    -0.541169741821289, -0.221188925170898, 2.2862017364502, 0.26431930847168,
    0.340397464707319, 0.279939501342773, -0.023313720703125,
    0.00999159458705358, 0.232991496276855, -0.0116082901000977,
    0.0931091720581054, -1.90108821258545, 0.02797646484375, -0.594937010192872,
    0.0314735229492188, -0.63121898803711, 0.567689099121094, 0.0116568603515626,
    -0.396566389160156, -0.104911743164062, 0.349997232055664, 0.832299829101562,
    0.338340371704101, 5.54174426651002, 10.1030008666992, 95.3158157226563,
    0.0525287269592284, 0.162613201904297, -0.155619085693359, -1.03818912506104,
    12.7153032714844, 2.01823965911865, -0.14129086151123, -0.209987410926819,
    -0.101123263549805, 0.116277182006836, 0.0828940812763412, 0.583571571350095,
    4.03604218597412, -4.88626443786621, 4.86484495697021, 0.544317094116211,
    -0.519875155857631, 1.09268494720459, -0.288762287521362, -2.64246453094482,
    0.0797123538746555, 0.282879215812683, 0.230441558074951, -2.85505652160645,
    0.109128783820657, -0.277027820354959, 0.345942671933382, 0.415752521584251,
    0.402598814392089, -4.94658869018555, -0.0490170977783206,
    -0.307713792514801, 0.0554949816022599, 0.875138790893555, 0.199758782512386,
    -0.415275650024414, -1.07078998911004, 0.973159272496841, 4.04871902160645,
    -7.75778627471924, 0.863481930541992, -8.34135784606934, -0.645104365808824,
    0.164444994245257, -0.853282177734375, 0.0710434956426204,
    -0.878386059134347, 0.361654092407226, -0.608238320486886,
    0.00779552536010719, 0.563084639282227, -0.215107084987475, 2.46236603851319,
    -0.303421217974494, 0.260530828857422, 0.12565871288593, -0.108346353803362,
    1.85475219268799, -0.220205377578735, -0.109154497483197, -0.376356777560327,
    5.7631517578125, -2.32204658203125, -0.353673626474234, -0.191803923034668,
    -0.439755056762696, -0.247685865431566, -0.207200692749023,
    -0.592213339937651, 0.175602274867467, 6.49840822448731, 0.0797037826538087,
    -1.54045409545898, -0.240557246985802, -1.52777725982666, 0.773942671966553,
    0.0807601856231687, 3.22035338287354, 0.991415972900395, -0.398462870671199,
    0.204453004237584, -0.589878765433176, -0.173596150016784, -8.25961411285401,
    -0.645972201919556, 0.442124773768375, -2.69407976355919,
    -0.0735727224496691, -0.458916020965576, 0.467731521606446,
    -0.134090321731567, -0.711651324462888, -0.280056069946291,
    -10.1148034378052, 0.186634660557338, -0.260418020531439, 4.28214765014648,
    0.463688048171997, -0.407144989929199, -0.344527474506085, 0.18222170628139,
    -0.083770437344638, -0.00971405029296888, 0.364043748779297,
    1.01782479100721, 3.46176961004084, 2.81868711730957, 1.45580248240595,
    0.47134514831543, -22.2369182281494, -13.0355755096436, 1.15571941955566,
    1.0577143661499, -0.885629965209944, 3.41446924987793, -1.23665966203962,
    3.30265346624201, 0.41539221862793, -3.2349660899571, 1.70947857055664,
    -1.03836397796631, 0.51416121270316, -1.89933968353271, 2.34302893066406,
    0.906903735351564, -0.261670022028143, 2.89486469970703, -1.73139346801758,
    -0.479321771327428, -2.0719195010376, -1.1546214185161, 0.611013763427735,
    -1.08933359985352, 0.0062072781372071, -25.2688676055908, -57.4925095184326,
    0.209998339233397, -4.52423149749756, 1.08285756632487, -3.75315036057692,
    -3.73310952758789, 3.80130216064453, -0.725726983337404, 1.14353800048826,
    -3.29007597885132, 23.2624305175781, 0.457813141290072, -1.12499110303606,
    28.0577714447021, 0.0646283238337588, 4.26990794677738, -1.19095810839844,
    -14.0232030029297, 1.10647378596256, -4.42168026855469, -0.162613201904293,
    -0.335997342773438, -1.87134573339844, 19.1706976812744, 1.29497700891495,
    3.29095024337769, -0.222063189697266, 0.34183742980957, 51.9584150665283,
    42.7113191711426, -186.299650717163, 66.0664217285156, -20.4088311035156,
    67.2519244262695, 1.58766437988281, -14.1823191467285, 2.22762601318359,
    17.9486507263184, 0.194960989379877, -15.9771842193603, 8.68756659851075,
    24.6385228820801, -27.9152663269043, -99.2806053497314, -0.192338195800781,
    -0.546415328979493, -0.0253536712646486, 0.421395501708984,
    -0.163487466430664, -0.373310952758789, 87.2402342926026, -51.9636606536865,
    8.1105520111084, -103.105512652588, 79.8369622833252, 3.98140065307617,
    60.3399890808105, 4.44738364562989, 3.03369790649414, -1.04824316711426,
    -3.56699926757812, 3.04856040344238, 41.3588319488525, -196.8581434021,
    15.9562018707275, 9.58019067993164, -2.09561206970215, 1.46351881713867,
    2.350023046875, 2.67612371520996, 4.44913217468262, 3.40263753662109,
    5.87855467529297, 4.28914176635742, 4.03997637634277, 1.44341073303223,
    12.859556918335, -0.396041830444336, 3.12374715270996, -46.7032109985352,
    2.78802957458496, 5.78063704833984, 1.04999169616699, -0.952074069213867,
    -24.4855265899658, 1.80273345336914, -29.0308278625488, -29.4801998291016,
    41.1507569915772, 15.9701901031494, 9.55046568603516, -3.28636035461426,
    -3.49443531188965, -2.77753840026856, -3.65967130737305, 30.2005937988281,
    -11.6399579040527, 4.20258957824707, -4.01462270507812, 1.72404964599609,
    31.9972074005127, -14.4087536590576, 12.6855782775879, 10.8767249725342,
    -23.4128040161133, 1.65410848388672, -5.75528337707519, 1536.40100511475,
    33.4091446105957, -49.8155927124023, 25.4953021179199, -9.92814796142578,
    -0.470354315185548, 1.86043491210938, 88.4948038879395, 20.8074957275391,
    -126.464986532593, -206.011692993164, -197.755138806152, -37.0828041503906,
    -71.362716229248, -171.89701690979, -10.0654074920654, -205.850828320313,
    -8.47686884765621, -219.336358639526, -69.0240586212158, -177.497555465698,
    -25.6500469390869, -189.119153814697, -80.33879012146, 197.906386569214,
    -267.380691421509, -156.331611282349, -51.381400479126, 235.431568569946,
    -291.755186416626, -7.12088456726077, -74.088673022461, -121.778928671265,
    -81.717505279541, -182.840185986328, -5.84096130065917, -253.797243475342,
    -9.01541579589842, -252.001504138184, 15.0924285186767, -241.870526806641,
    -131.419443603516, -36.8248961151123, -249.995941314697, 41.343095187378,
    -150.05963757019, -123.484618762207, 26.1055387573242, -200.622726452637,
    -71.5655455993652, -132.863728601074, 9.42107453613279, -110.973893389893,
    -77.7002597808838, -243.205528738403, -21.0706493499756, -90.7399151916504,
    -115.665196838379, -102.836239178467, -236.598711712647, -13.7932714324951,
    -77.1232451934814, -225.518283105469, -72.6872269866944, -38.857561138916,
    5.82872159729003, -115.202710903931, -59.9710494506836, -112.653355545044,
    -128.480166265869, -47.1849307525635, -179.689336633301, -63.1918399658203,
    -69.3693931091309, -282.126036923218, -51.6856445343017, -103.09065015564,
    -44.736115814209, -108.437651998901, -135.803005938721, -153.627511102295,
    -60.9537227783203, -61.4878984039307, -157.316907403564, -135.129822253418,
    -33.7929467376709, -95.9269266265868, -193.508836001587, -169.145706445313,
    7.01772135314944, -72.0332771209716, -187.811254083252, -196.093161941528,
    33.9188408294678, -141.327483480835, -131.276064221191, -190.070353619385,
    81.9421912628174, -151.19792998352, -140.001224194336, -225.374903723145,
    12.8997730865478, -112.299278411865, -187.157304217529, -164.825091156006,
    -116.645247372437, 67.2868950073243, -229.964792486572, 408.710797695923,
    149.387328149414, 17.1355847167969, -37.5155650909424, -52.2329341278077,
    -97.1255432922363, 69.9630187225342, -86.7690057128904, -88.1704517486571,
    -12.6392422576905, 206.269601028442, 367.692054647827, -135.907043417358,
    458.807029321289, -202.966629647827, 13.2879465362549, -313.453557696533,
    -150.136572848511, -21.4798051483154, -53.421933883667, -26.5723960144043,
    175.937867550659, -60.0969435424805, -95.0491650421144, -139.570211782837,
    518.192321539307, -126.506076965332, -89.1024177337647, -114.606462496948,
    259.173096047973, -174.383425222778, 338.362228317261, -171.021878118897,
    -66.0428165863037, -161.152305880737, -120.388848074341, 89.6794323211669,
    34.0963165283204, 170.574254681397, 241.810202554321, 292.465963476563,
    -53.145666293335, 342.653992877197, -165.84448359375, -133.782580618286,
    77.9188259124756, 26.6738106994628, -212.576545321655, -237.380304199219,
    38.7412839569092, 165.343530020142, -121.347041995239, -66.3951451904298,
    792.797935006714, -391.992237158203, 921.505410049438, 0.0268746731853623,
    0.0298055364674886, 10.5331512918158, 40.3563553484106, 10.5196057591678,
    20.2662766634931, 6.45840157690791, 6.92889450998839, 1.12571940116754,
    2.10476284038953, 3.08719693673889, 2.34711734374375, 6.712959053828,
    1.67516356137947, 3.52578266489849, 2.7950869454372, 5.30635133801325,
    8.09240434150817, 0.0246115914116437, 0.0242304215956303, 0.0286001704913705,
    0.0291879002613391, 0.0290566979699974, 0.0256623134360342, 5.83900879046717,
    8.41617526473238, 5.39441247987161, 7.87686072741409, 6.68963677370888,
    2.52977274804707, 10.6956000922393, 4.41072504738288, 0.428610002607764,
    0.463463452687495, 0.267164360939963, 0.133570538898035, 5.47880969017152,
    19.9030278401895, 29.5682839429925, 2.38941495671309, 0.466781318594176,
    0.148329886727391, 0.184086280646204, 0.123922749514189, 0.140527423348506,
    0.172344359457348, 0.203513771720173, 0.146341971213772, 0.110380285708005,
    0.180868345187154, 1.08877056354086, 0.21336181207467, 0.148820077386298,
    2.7841223440841, 1.80721584603855, 0.122362435454045, 0.277264245692295,
    0.572730719684492, 2.83846958304864, 0.791002224241396, 3.41240518071897,
    2.57254123743482, 3.39147881639427, 2.43661905639406, 1.14417187856408,
    0.388130879016578, 0.50287990019218, 0.207715646077542, 0.248161738823364,
    2.99815178411407, 0.377406159963217, 1.19324898120731, 0.247367734107045,
    0.696800214597851, 3.3574195521139, 2.31088967942236, 3.06883201574632,
    0.917886646062514, 1.25601183434853, 3.57924996802984, 1.21027032387873,
    69.7311566807075, 27.2723306188125, 26.3303825660204, 3.08578186879948,
    1.38591370799485, 2.2362401297356, 1.54618941760835, 11.2602237974904,
    27.7424962272088, 25.7357686226883, 27.6881360688327, 25.7886113936362,
    26.9015625573976, 26.1576633506573, 21.3623414912971, 22.3664704290412,
    20.7971996349278, 23.6078990430166, 20.6784744922384, 24.8006568736018,
    22.7942708627011, 23.3636201100565, 21.9074784705648, 38.3531917663255,
    43.0712757026176, 25.6554588377047, 31.5935230157821, 45.4776019124526,
    42.8369261738649, 25.4369306825786, 26.0756476905114, 27.5281213825806,
    29.0962862516123, 27.4248549999529, 25.3614468595669, 27.0865080499951,
    24.6604111795759, 25.8864725793388, 27.3497285975656, 27.3079817484561,
    26.6810983765145, 27.1678926426133, 25.6042308562437, 25.9581905759086,
    25.9891178389962, 25.9788539470771, 27.3201717528639, 24.7470892919295,
    24.0885107164743, 26.4817272846584, 26.8516637316326, 26.3920518955321,
    24.9853480560993, 26.5529845590576, 28.6171543292616, 29.8570802836042,
    29.6968437256118, 26.6325071550816, 27.7073314750816, 28.2453797754439,
    31.1678392491771, 26.9398833099661, 29.7013822479879, 29.5148316527414,
    27.9715395012743, 26.0075759783692, 27.6823097699155, 26.5883490230618,
    28.760241403099, 27.5392667641776, 26.9859966550091, 26.8965842040549,
    26.8769353575258, 26.2562818922336, 29.1601250757688, 29.820672511046,
    28.3921898532625, 24.9615515405967, 26.147938740201, 28.1061337504931,
    24.8047004984786, 27.999245781427, 26.9655385644042, 24.7520411913248,
    32.0013014236074, 30.1638468304317, 30.2766202441249, 28.3901162807007,
    27.5831330789031, 28.2199694101508, 28.4790909087028, 28.9602594586738,
    29.0282784094054, 26.7968383251726, 28.1241859021578, 26.5424676967736,
    28.2531453827206, 26.8843210748958, 30.2993535214809, 28.3033550922335,
    28.3060869431279, 27.4115739339638, 30.0336882554787, 28.7776257726441,
    27.4625123568408, 26.9074007401214, 26.7849031477795, 65.9778067654654,
    73.3821836406903, 73.6231496726087, 78.5679413749248, 71.5014802006349,
    74.2468241072797, 72.0857761769689, 69.1172527575299, 71.9709065997591,
    70.9527197097668, 72.5818522297698, 75.1280369358456, 62.5997132415406,
    71.0045551289196, 73.3813624756525, 73.0761115273891, 76.5037491235949,
    71.2411868705645, 74.6675547196893, 72.8395662578339, 71.1310921361007,
    74.3476770521163, 75.5608426560147, 73.586454321016, 77.2358443450363,
    73.8399674046365, 78.0913200663529, 74.245337592819, 75.0896349543626,
    73.2758082092055, 76.7646455350516, 76.2313179495815, 77.4352700079064,
    77.5205229926473, 77.5558136199105, 80.4106430994672, 80.416896705231,
    80.8889525069117, 73.0977565528943, 77.7823183151069, 78.2524375374584,
    75.1040357926508, 70.6210538334645, 72.6639648034683, 76.2880185959108,
    74.7817190216895, 76.2637688563543, 76.2351385724043, 74.5196207340348,
    73.8506087684824, 73.1643438361246, 75.4674734076707, 77.7225164450469,
    76.4677096924535, 90.3855571560722, 84.7177601022202, 99.773373030477,
    0.0472102844238281, 0.0708154266357422, 50.2142573364258, 162.469822521973,
    22.5874983032227, 79.0536212677002, 10.2673625976563, 21.3268088562012,
    2.28795026550293, 4.20084104919434, 5.92751348876953, 5.78675690002441,
    14.3274470581055, 3.24439565734863, 10.8268918945313, 11.0603205230713,
    14.8869763549805, 12.5544385986328, 0.0472102844238281, 0.0472102844238281,
    0.0515816070556641, 0.0751867492675781, 0.0515816070556641,
    0.0515816070556641, 36.9647784393311, 35.1917699798584, 24.5510964294434,
    22.14774324646, 21.6581551116943, 7.56501094665527, 38.8785434875488,
    11.904860055542, 0.938085836791992, 1.98807753295898, 0.401287417602539,
    0.401287417602539, 34.164509161377, 75.2086058807373, 52.7845950439453,
    13.0011877716064, 0.681926330566406, 0.401287417602539, 0.354951397705078,
    0.331346255493164, 0.401287417602539, 0.378556539916992, 0.799077777099609,
    0.424892559814453, 0.402161682128906, 0.331346255493164, 5.74479220275879,
    0.542044006347656, 0.401287417602539, 5.6512458984375, 5.51136357421875,
    0.424892559814453, 0.542044006347656, 2.35352010498047, 7.37004995727539,
    3.35630151672363, 6.22651195678711, 3.84676391601563, 16.166899621582,
    9.5627053894043, 4.26641088867188, 0.673183685302734, 1.09370492248535,
    0.440629321289063, 0.486965341186523, 15.4902188781738, 0.790335131835938,
    5.17652026062012, 0.417024179077148, 2.60967961120605, 8.11667186279297,
    9.09584813232422, 13.5301178100586, 3.14647803039551, 2.05015031433106,
    13.0169245330811, 2.56334359130859, 310.305331137085, 96.2469074432373,
    45.3332384857178, 12.1059408966064, 4.82681445007324, 8.95684007263184,
    5.27006656494141, 35.6726154693604, 55.0026041473389, 50.5928138763428,
    48.2130658355713, 46.6726117401123, 46.5318551513672, 49.7054353820801,
    47.8388806182861, 38.7159302856445, 40.5588799072266, 38.202737008667,
    37.2690224945068, 46.9987124084473, 45.6226200439453, 39.7886528594971,
    41.9358465362549, 90.7949938568115, 86.3624727081299, 35.5423500549316,
    87.855716519165, 93.5253219726562, 86.9228762695312, 39.812258001709,
    56.6593354248047, 57.7093271209717, 49.6127633422852, 50.779032220459,
    48.4919562194824, 46.6254014556885, 44.3628048614502, 55.0725453094482,
    53.9526124511719, 55.8182929504395, 47.2094101593018, 50.7326962005615,
    51.7118724700928, 54.1388307952881, 45.9959309967041, 42.3318883666992,
    55.7247466461182, 48.9588134765625, 45.4591325775147, 63.0056216217041,
    53.5556963562012, 62.8657392974854, 44.4327460235596, 47.7226034362793,
    65.3853696624756, 54.7228394989014, 55.5857385864258, 55.2596379180908,
    54.6284189300537, 53.6955786804199, 52.5284355377197, 46.6490065979004,
    47.1858050170898, 59.5994870269775, 46.8124940643311, 60.9056382293701,
    57.1261926818848, 45.3883171508789, 64.079218460083, 48.188586428833,
    47.6063262542725, 47.3492924835205, 53.5329654785156, 49.8925279907227,
    66.9022186157227, 63.5660251831055, 59.0854194854736, 39.9984763458252,
    54.6520240722656, 56.3550913696289, 39.7659219818115, 56.8455537689209,
    57.312411026001, 53.3222677276611, 69.3991181030274, 66.5525128051758,
    66.7396054138184, 52.6927972686768, 55.8191672149658, 61.4424366485596,
    63.9620670135498, 53.6020323760986, 49.0060237609863, 47.0459226928711,
    52.5756458221436, 41.1892246307373, 56.0753267211914, 52.879015612793,
    57.9663608917236, 67.8123279876709, 59.879251675415, 55.3986459777832,
    65.7586806152344, 61.4660417907715, 49.2621832672119, 51.6655364501953,
    56.8927640533447, 116.125059979248, 143.378508059692, 129.774952029419,
    169.978006274414, 116.614648114014, 137.895120950317, 121.981758041382,
    127.347993704224, 145.432155432129, 151.941054830933, 132.038422888184,
    200.521311767578, 140.322079275513, 122.751985089111, 175.018141268921,
    155.441609994507, 130.218204144287, 113.97786630249, 125.038186825562,
    122.868262271118, 123.381455548096, 149.958222885132, 160.784240515137,
    133.157481481934, 149.514970770264, 132.551616165161, 154.85847555542,
    135.117582550049, 132.831380813599, 135.024910510254, 157.635139691162,
    143.868096194458, 151.521407858276, 141.442012133789, 144.171465985107,
    155.651433480835, 157.587929406738, 155.511551156616, 151.497802716064,
    157.495257366943, 164.611770611572, 207.777707336426, 142.514734707642,
    151.218038067627, 140.461961599731, 129.587859420776, 139.387490496826,
    152.268029763794, 136.821524111939, 131.711447955322, 151.33518951416,
    154.601441784668, 158.008450643921, 161.741560171509, 202.178043045044,
    195.784546563721, 269.167688113403, -0.0463360198974609, -0.0690668975830078,
    -34.9976832550049, -88.3147053955078, -36.9586585876465, -39.6881124389648,
    -31.7314309844971, -19.4585055633545, -4.43164688415527, -10.2420089263916,
    -7.76871458129883, -3.59147867431641, -18.8045556976318, -6.15831932373047,
    -14.9315638458252, -8.1184203918457, -12.9015216156006, -34.7651288909912,
    -0.0463360198974609, -0.0690668975830078, -0.0646955749511719,
    -0.041964697265625, -0.0646955749511719, -0.041964697265625,
    -0.0646955749511719, -31.540841317749, -12.3376209960938, -21.8802183013916,
    -12.7109319488525, -4.3109983795166, -21.2507478424072, -11.7308814147949,
    -0.973930682373047, -0.881258642578125, -1.1138130065918, -0.158241879272461,
    -4.94046883850098, -45.5867751983643, -100.327099987793, -3.47083016967773,
    -2.02392237854004, -0.438006527709961, -0.438006527709961,
    -0.274519061279297, -0.180972756958008, -0.717771176147461,
    -0.204577899169922, -0.274519061279297, -0.157367614746094,
    -0.577888851928711, -1.04387184448242, -0.554283709716797,
    -0.204577899169922, -12.7808731109619, -4.84692253417969, -0.181847021484375,
    -0.866396145629883, -1.49586660461426, -12.0421195861816, -2.12533706359863,
    -13.2092627288818, -12.2056070526123, -2.82562294921875, -5.2517070098877,
    -2.84922809143066, -1.02900934753418, -1.51947174682617, -0.586631497192383,
    -0.656572659301758, -2.49952228088379, -1.37871515808105, -2.05539590148926,
    -0.820060125732422, -2.40510171203613, -7.04919487609863, -6.27896782836914,
    -5.34612757873535, -1.12255565185547, -4.90200119934082, -7.46884184875488,
    -5.04188352355957, -2.96637953796387, -68.8588226257324, -128.61480300293,
    -5.13630409240723, -4.66944683532715, -4.20258957824707, -3.12899273986816,
    -47.6946269714356, -48.1763467254639, -46.3097919616699, -50.6033050506592,
    -43.2996991973877, -49.5987751098633, -45.8429347045898, -39.5893205474854,
    -42.1797663391113, -49.5533133544922, -52.3527083679199, -46.052758190918,
    -46.052758190918, -47.569607144165, -44.1162622650147, -46.9165315429688,
    -60.636364755249, -57.6726080108643, -64.2523228363037, -63.8090707214355,
    -56.5063391326904, -49.4597670501709, -52.4925906921387, -49.5297082122803,
    -45.7030523803711, -53.7996161590576, -47.9892541168213, -41.9926737304688,
    -45.1662539611816, -47.4297248199463, -53.310028024292, -50.0665066314697,
    -42.4595309875488, -53.4726412261963, -47.1499601715088, -40.1497241088867,
    -49.0628509552002, -48.7830863067627, -52.0729437194824, -47.336178515625,
    -40.8727408721924, -45.2598002655029, -42.7165647583008, -44.7929430084229,
    -37.7227657836914, -47.3134476379395, -48.7830863067627, -47.2199013336182,
    -51.536145300293, -52.6560781585693, -49.6223802520752, -55.0594313415527,
    -56.2956413818359, -57.1830198760986, -53.7760110168457, -60.7762470794678,
    -53.9631036254883, -53.6833389770508, -54.8496078552246, -51.1628343475342,
    -52.7959604827881, -50.486153604126, -47.4533299621582, -51.769573928833,
    -45.6794472381592, -51.6769018890381, -49.6223802520752, -70.553147277832,
    -51.8631202331543, -55.1066416259766, -48.2226827453613, -49.4597670501709,
    -49.6468596588135, -48.0364644012451, -52.6097421386719, -45.002766494751,
    -41.9935479949951, -62.5728606811523, -45.9364810089111, -56.9495912475586,
    -55.7597172271729, -53.1928765777588, -56.7397677612305, -52.4925906921387,
    -52.3054980834961, -59.0259694976807, -52.0257334350586, -55.4563474365234,
    -50.3698764221191, -52.8895067871094, -46.286186819458, -54.7097255310059,
    -53.4962463684082, -58.7462048492432, -52.5398009765625, -52.7496244628906,
    -55.5026834564209, -55.3628011322022, -44.6294555419922, -49.9493551849365,
    -119.726155563355, -121.91968526001, -108.876532791138, -120.21574369812,
    -114.47619708252, -107.313347817993, -122.269391070557, -122.456483679199,
    -112.516970278931, -131.089845877075, -115.57339906311, -111.000121325684,
    -105.983591473389, -113.706844299316, -134.846560546875, -137.436132073975,
    -118.886861618042, -127.893534768677, -112.936617251587, -114.010214089966,
    -113.800390603638, -114.87311317749, -124.556467071533, -116.855945123291,
    -120.542718630981, -120.286559124756, -127.776383322144, -114.849508035278,
    -118.62982784729, -112.633247460938, -116.202869522095, -130.296013687134,
    -121.966895544434, -129.269627133179, -118.420004360962, -126.749996768188,
    -138.089207675171, -134.729409100342, -123.926122348022, -131.836467782593,
    -114.87311317749, -125.022450064087, -108.550432122803, -118.769710171509,
    -137.92572020874, -140.795930648804, -138.789493560791, -143.7832925354,
    -111.909356433105, -135.266207519531, -111.956566717529, -126.61011444397,
    -127.099702578735, -123.063223260498, -145.486359832764, -148.799822387695,
    -148.729881225586, -0.222937454223633, 0.305992584228516, 26.8626518371582,
    56.0215594528198, -187.875949658203, 55.0524372253418, -20.6785417098999,
    66.4226845230103, 1.83158418273926, -12.3170757797241, 2.69317187347412,
    18.5776840530395, 0.975679211425778, -15.8032055786133, 9.44467967834473,
    26.5391739624023, -28.6509599258423, -97.986693850708, -0.205015031433105,
    -0.525870112609863, -0.0305992584228516, 0.404347343444824,
    -0.15693048248291, -0.331346255493164, 86.2431356002808, -50.4874650009155,
    10.2280206939697, -101.781001895142, 80.1114813446045, 1.61957503509521,
    63.2740208312988, 8.88689891052246, 2.46892302246094, -1.12299278411865,
    -3.87517751312256, 2.78715531005859, 39.4166533035278, -205.659364389038,
    14.8187837219238, 10.6131342178345, -2.43919802856445, 1.47051293334961,
    2.36838260192871, 2.67087812805176, 4.37394542541504, 3.25750962524414,
    5.52316614532471, 4.16761899719238, 4.01112564697266, 1.31008539276123,
    12.5046055206299, -0.57657745513916, 3.10713612670898, -48.4578599029541,
    1.77956544342041, 5.62370656585693, 1.13086116485596, -1.17064020080566,
    -24.3071766265869, 2.63372188568115, -29.0045999267578, -29.6401902374268,
    33.8156776153565, 12.0298798828125, 9.72881564941406, -3.60983822937012,
    -2.96637953796387, -2.57602042694092, -3.59803565826416, 30.9035024780273,
    -11.2867550354004, 4.38093954162598, -3.87124332275391, 1.58722724761963,
    28.5818930282593, -19.3632107299805, 12.6886382034302, 10.7748731552124,
    -23.4560801101685, 0.864210484313968, -6.26541672821045, 1462.47275963287,
    -14.1062581329345, -40.5142924163818, 22.0104837158203, -11.8969916748047,
    -0.000874264526365964, 1.71224707489014, 86.6549141921997, 32.8758432495117,
    -138.803481793213, -221.849869152832, -192.884611129761, -22.5993008743286,
    -59.3297764205933, -183.802751229858, -39.8328032180786, -213.883570788574,
    -37.70746615448, -219.214398738098, -54.3775050109863, -180.980625338745,
    -50.4240808227539, -207.313910005188, -75.2352709487914, 176.736945327759,
    -285.621783631897, -173.896022749329, -55.5293485244751, 235.483587309265,
    -277.28348571167, -13.9759927185059, -73.7682550735474, -94.2990460784912,
    -99.2347064620972, -186.836886268616, -0.539421212768552, -271.863045648193,
    -28.585827218628, -262.624255265808, -15.7712949234009, -227.784376757812,
    -137.445748983765, -12.8219635437012, -234.964711312866, 24.3390872817993,
    -146.332647894287, -128.43776443634, 43.3770716079712, -167.472801274109,
    -63.0091186798096, -133.686848652649, -5.61190399475098, -112.870173147583,
    -64.9552315155029, -250.130578051758, -21.8937694015503, -123.902517205811,
    -88.4305454452515, -112.957162467957, -231.355310215759, 22.4244479690552,
    -46.1904548538208, -215.946397938538, -90.426928491211, -36.157395149231,
    -20.5439049728394, -113.003935620117, -40.7398526641846, -142.373978118897,
    -116.400453305054, -28.1985280334473, -204.615492544556, -71.7736205566406,
    -56.2628564620972, -270.710765002441, -73.1116824142456, -75.6693432861328,
    -73.2664272354126, -103.742414360046, -90.531840234375, -165.697170021057,
    -68.4033308074951, -28.7117213104248, -160.531578067017, -137.539295288086,
    -29.3661083084106, -104.986492781067, -204.354087451172, -166.865187428284,
    17.6863713684082, -53.6649794219971, -179.406512059021, -175.613078279114,
    35.1139604370117, -110.756201522827, -135.84409637146, -201.615016690063,
    79.4859450759888, -112.402004493713, -147.614319689941, -214.275241296387,
    18.4465443740845, -72.9897225128174, -199.262370849609, -155.254080253601,
    -93.6923064971924, 55.6023496124268, -253.527095736694, 429.854885266113,
    190.573929986572, 58.0306193344117, 3.26319234466553, -41.6587046813965,
    -45.2287638748169, 80.0585883407591, -83.9416342346191, -141.937720120239,
    -72.8231751205445, 227.343747436523, 394.24915103302, -102.466425283813,
    400.174915992737, -191.949585218811, -108.436777734375, -350.36282033844,
    -193.999298400879, 35.6791724533081, 29.7879409423828, 32.0212496749878,
    204.268409527588, -63.4532450592041, -84.3354904037476, -119.522889060974,
    484.071962736511, -101.932249658203, -25.5272127731323, -75.6820201217651,
    242.84402035675, -177.914142512512, 354.267722845459, -209.960745858765,
    -70.2152440383911, -59.8477781524658, -86.6807049957275, 1.08758507080069,
    -62.5920945007325, 224.349391433716, 219.811521409607, 273.66403057251,
    -54.0977403625488, 361.849781950378, -72.9748600158694, -97.6243112045288,
    -2.5528524169921, -15.8233136627197, -185.272827030945, -209.749610975647,
    95.5147109024049, 136.601209451294, -127.141230143738, -54.9230460754395,
    767.240122975159, -306.530694044495, 959.133755264282, 0.0, 0.0, 0.0, 2.0,
    0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 3.0, 0.0, 1.0, 3.0,
    5.0, 0.0, 2.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 1.0, 0.0, 0.0, 1.0,
    1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 2.0, 4.0, 3.0, 1.0,
    1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 1.0, 1.0, 2.0, 0.0, 1.0, 1.0, 0.0, 1.0, 2.0, 1.0, 2.0, 2.0, 3.0, 1.0,
    1.0, 2.0, 1.0, 1.0, 0.0, 0.0, 2.0, 0.0, 2.0, 1.0, 1.0, 3.0, 2.0, 2.0, 2.0,
    0.0, 0.0, 1.0, 2.0, 8.0, 9.0, 9.0, 9.0, 9.0, 8.0, 9.0, 9.0, 9.0, 7.0, 9.0,
    9.0, 8.0, 8.0, 8.0, 7.0, 9.0, 9.0, 10.0, 9.0, 10.0, 11.0, 8.0, 9.0, 9.0, 8.0,
    9.0, 9.0, 9.0, 8.0, 8.0, 9.0, 10.0, 7.0, 10.0, 8.0, 10.0, 9.0, 9.0, 9.0, 9.0,
    8.0, 8.0, 9.0, 6.0, 9.0, 9.0, 9.0, 10.0, 9.0, 8.0, 10.0, 8.0, 8.0, 10.0, 8.0,
    8.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    18.0, 39.0, 0.0, 0.0, 37.0, 42.0, 0.0, 3.0, 21.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 17.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 20.0, 37.0, 22.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 19.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    19.0, 0.0, 18.0, 18.0, 36.0, 0.0, 0.0, 18.0, 0.0, 0.0, 0.0, 0.0, 18.0, 0.0,
    18.0, 0.0, 0.0, 36.0, 18.0, 19.0, 37.0, 0.0, 0.0, 0.0, 3.0, 39.0, 42.0, 44.0,
    39.0, 39.0, 40.0, 41.0, 40.0, 45.0, 35.0, 38.0, 37.0, 41.0, 41.0, 46.0, 38.0,
    38.0, 42.0, 39.0, 41.0, 40.0, 44.0, 40.0, 40.0, 40.0, 40.0, 39.0, 42.0, 43.0,
    40.0, 39.0, 39.0, 40.0, 40.0, 39.0, 38.0, 38.0, 38.0, 46.0, 38.0, 38.0, 39.0,
    39.0, 40.0, 38.0, 42.0, 38.0, 45.0, 38.0, 37.0, 39.0, 47.0, 38.0, 38.0, 43.0,
    35.0, 42.0, 0.0, 0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 9.0, 13.0, 0.0, 0.0, 12.3333333333333, 8.4, 0.0, 1.5, 10.5, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 8.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 10.0, 9.25, 7.33333333333333, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    9.5, 0.0, 0.0, 0.0, 0.0, 0.0, 9.5, 0.0, 9.0, 9.0, 12.0, 0.0, 0.0, 9.0, 0.0,
    0.0, 0.0, 0.0, 9.0, 0.0, 9.0, 0.0, 0.0, 12.0, 9.0, 9.5, 18.5, 0.0, 0.0, 0.0,
    1.5, 4.875, 4.66666666666667, 4.88888888888889, 4.33333333333333,
    4.33333333333333, 5.0, 4.55555555555556, 4.44444444444445, 5.0, 5.0,
    4.22222222222222, 4.11111111111111, 5.125, 5.125, 5.75, 5.42857142857143,
    4.22222222222222, 4.66666666666667, 3.9, 4.55555555555556, 4.0, 4.0, 5.0,
    4.44444444444445, 4.44444444444445, 5.0, 4.33333333333333, 4.66666666666667,
    4.77777777777778, 5.0, 4.875, 4.33333333333333, 4.0, 5.71428571428571, 3.9,
    4.75, 3.8, 4.22222222222222, 5.11111111111111, 4.22222222222222,
    4.22222222222222, 4.875, 4.875, 4.44444444444445, 6.33333333333333,
    4.66666666666667, 4.22222222222222, 5.0, 3.8, 4.11111111111111, 4.875, 4.7,
    4.75, 4.75, 4.3, 4.375, 5.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 1.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 3.0, 2.0, 0.0, 0.0, 0.0, 2.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 5.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 4.0, 6.0, 9.0, 8.0, 6.0, 8.0, 8.0, 4.0, 7.0, 6.0, 8.0,
    4.0, 4.0, 4.0, 7.0, 7.0, 7.0, 4.0, 8.0, 8.0, 7.0, 7.0, 6.0, 6.0, 7.0, 5.0,
    8.0, 8.0, 7.0, 8.0, 8.0, 7.0, 10.0, 9.0, 8.0, 8.0, 8.0, 9.0, 5.0, 8.0, 10.0,
    6.0, 7.0, 8.0, 6.0, 4.0, 9.0, 7.0, 8.0, 6.0, 5.0, 7.0, 9.0, 9.0, 8.0, 5.0,
    10.0, 50.0, 50.0, 50.0, 48.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 43.0, 49.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 48.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 46.0, 44.0, 41.0, 42.0, 44.0, 42.0,
    42.0, 46.0, 43.0, 44.0, 42.0, 44.0, 46.0, 46.0, 42.0, 43.0, 43.0, 46.0, 42.0,
    42.0, 43.0, 43.0, 44.0, 44.0, 43.0, 45.0, 42.0, 42.0, 43.0, 42.0, 42.0, 43.0,
    40.0, 41.0, 42.0, 42.0, 42.0, 41.0, 45.0, 42.0, 40.0, 42.0, 43.0, 42.0, 44.0,
    46.0, 41.0, 43.0, 42.0, 44.0, 45.0, 43.0, 41.0, 41.0, 38.0, 42.0, 38.0,
    0.00759315567054217, 0.0101712848314297, 81.869803757427, 87.7856289824233,
    7.08372162339915, 12.9570641125112, 3.58763054707985, 18.3560886158295,
    0.582607067008958, 3.41800043537413, 2.33856911623612, 0.453448369698706,
    4.88651157353227, 8.05576862942919, 11.1527116852074, 3.79428022390292,
    3.14786353688984, 14.6262226329465, 0.00511418216419807, 0.0205736363516741,
    0.00326790058635369, 0.00551790426158598, 0.0466105588774311,
    0.0673334081698928, 20.8953531794399, 63.7904805598862, 0.316083649177074,
    1.20077707245832, 1.88840628569388, 1.14520676787514, 32.6659148815703,
    2.9762866308599, 0.0503823172525437, 0.907214812901224, 0.300905233603951,
    0.112506531860787, 20.8541721685658, 37.8619727407378, 157.267508374022,
    0.17642633234815, 1.38580711487717, 0.267601013173333, 0.0872212231012922,
    0.321858555119944, 0.0587298957965909, 0.0126984331039997, 0.278239086713511,
    0.233155594696771, 0.190885085868667, 1.706870040868, 0.554155818555285,
    0.171868360660806, 0.419603630516237, 3.43623483183046, 1.66423061349525,
    0.160826448139235, 1.1440931539575, 0.511851782329893, 3.80998566727015,
    6.60138235039343, 2.0389022575149, 2.23073891907893, 6.12478078730715,
    2.54817690382473, 2.01614531686675, 0.0904157352780601, 0.690767246317523,
    0.417379895514757, 0.286704781742551, 0.820664862079022, 0.945281562705025,
    2.15115469969747, 0.438487386557205, 6.45843550088262, 1.65569884325194,
    4.68396246078608, 9.82558656640872, 3.54036099112595, 1.90460791817505,
    11.3737902025139, 5.57029004337448, 28.2828236921565, 18.731278888487,
    37.4252185499817, 11.9289436129155, 2.58804038877349, 1.93164204551492,
    2.16018398108628, 6.19705425366311, 11.9334711616036, 14.9548487996808,
    16.3238124816371, 12.4481822446272, 0.238157973572442, 17.7398931914487,
    7.6514294754376, 10.4338858772616, 0.0792108074448024, 13.4959744510201,
    36.1969348503363, 3.63554529294312, 2.37081764149915, 3.13599137588508,
    18.5632066320313, 4.21253880858347, 16.2878656147344, 8.86205306729983,
    28.8024799845014, 15.1925684919646, 7.78208960803073, 23.0027658149171,
    9.33414281881565, 23.4049123145219, 32.9901928549676, 14.0684178695519,
    15.0387231259048, 13.4085135755638, 29.8395137569748, 19.713638420806,
    8.47164182933756, 19.1066793326921, 46.6653736078643, 0.643284588794671,
    16.0021293082393, 12.3481209351066, 17.0784363527846, 10.5768680419127,
    6.14076743540597, 8.78070330719093, 12.6646744768697, 19.6444597559453,
    11.8532691470504, 7.67208845229992, 8.89618783311244, 17.883257631707,
    13.6596662841066, 0.461034270856885, 20.8966088374926, 21.2553218278806,
    16.2012490538483, 20.1819351781945, 41.2579066138441, 8.55468843530995,
    28.4268931248783, 1.27422084751621, 16.9503352333147, 16.2922927854822,
    5.30503862828933, 1.66647983858534, 1.76093172625745, 11.9818989977435,
    24.5088290794041, 19.1797626762838, 23.9979287371437, 9.09895141552063,
    7.44445568674504, 27.2090727891849, 4.69599686553243, 22.9225561620122,
    17.2000078950149, 10.7652535547059, 14.5548609584167, 7.71004107609423,
    3.89819682255758, 8.50728349572391, 13.6556388534631, 2.37795366837254,
    2.29235546687926, 4.24815837973889, 3.21358370206623, 3.39108600086736,
    7.7193228935991, 16.7550664034993, 0.664047644557705, 5.07181024761844,
    12.8085869080566, 6.43427219072122, 3.59253359562911, 2.34926003052021,
    5.91268217350811, 11.0747167932707, 1.02555751195024, 7.68196497801228,
    21.2782599827717, 4.52910134690202, 5.71661305752306, 9.56007753655346,
    26.5080348754515, 19.733622478238, 39.9549345546813, 30.1254400597832,
    3.44133882713186, 61.1307866050661, 33.6794191193304, 25.9738718016959,
    75.0380867867652, 78.3084053671662, 1.98301743557516, 25.9505322031018,
    19.2753146811559, 24.1648973759232, 20.7132958432096, 6.40487320179048,
    54.8736150747965, 20.2053492203019, 73.0647511849612, 57.5520152663025,
    42.0633255270404, 36.9236305801397, 24.9973838261109, 25.7633356998023,
    31.0556953284495, 38.994311103535, 2.02300868503141, 13.4151212631644,
    27.4013050959018, 43.5525828877535, 30.920841077884, 31.2152518201468,
    91.2923748868256, 2.04706879908119, 26.6612202747941, 83.3817476215429,
    82.5203173480609, 3.65040538061649, 65.927083911124, 26.2222495778569,
    43.1014912115068, 30.2934687304324, 3.55614029141267, 14.9145956767882,
    8.54994505319164, 83.5067212359636, 34.1500760156331, 79.0241156836383,
    70.8000649505987, 0.630428812915355, 38.1861498961801, 42.5818268722957,
    39.9703559679356, 54.5380823494675, 62.3043751516257, 52.1795652315624,
    85.6839036318491, 104.655010354868, 0.000361578841454389,
    0.00040685139325719, 2.3391372502122, 21.9464072456058, 0.236124054113305,
    -12.9570641125112, -1.19587684902662, 18.3560886158295, -0.291303533504479,
    -3.41800043537413, -0.101676918097223, 0.0266734335116886, 1.62883719117742,
    -8.05576862942919, -5.57635584260369, -3.79428022390292, -0.449694790984263,
    14.6262226329465, 0.00255709108209903, 0.00121021390303965,
    0.000217860039090246, -0.000367860284105732, -0.00517895098638123,
    -0.0336667040849464, 6.96511772647996, -63.7904805598862, 0.0117068018213731,
    -0.120077707245832, -0.209822920632654, -0.572603383937572,
    -10.8886382938568, -0.0783233323910501, 0.0100764634505087,
    -0.907214812901224, -0.300905233603951, -0.112506531860787, 1.04270860842829,
    -37.8619727407378, -52.4225027913407, -0.088213166174075, -0.461935704959056,
    0.014866722954074, -0.0872212231012922, -0.321858555119944,
    -0.0039153263864394, -0.00634921655199983, -0.039748440959073,
    0.116577797348385, 0.00867659481221213, -1.706870040868, 0.554155818555285,
    -0.171868360660806, -0.0466226256129152, -0.429529353978808,
    -0.128017739499635, -0.160826448139235, -0.381364384652501,
    -0.102370356465979, -0.76199713345403, 0.150031417054396,
    -0.0926773753415862, -1.11536945953946, 0.40831871915381, -0.637044225956183,
    -2.01614531686675, -0.0904157352780601, -0.230255748772508,
    0.417379895514757, -0.286704781742551, -0.205166215519756,
    -0.945281562705025, -0.717051566565824, -0.438487386557205, 3.22921775044131,
    0.413924710812986, 0.312264164052405, -4.91279328320436, -1.77018049556298,
    -0.952303959087527, 0.516990463750631, -5.57029004337448, -9.42760789738551,
    9.36563944424351, -9.35630463749541, 0.397631453763849, -2.58804038877349,
    -1.93164204551492, -0.0939210426559251, -3.09852712683156, 1.08486101469124,
    0.498494959989362, 5.44127082721237, 1.13165293132974, -0.0070046462815424,
    0.886994659572437, 3.8257147377188, -3.47796195908722, -0.0066009006204002,
    -0.421749201594377, 2.58549534645259, -0.106927802733621, 0.107764438249961,
    0.627198275177015, -1.54693388600261, 0.191479036753794, 3.25757312294688,
    -1.10775663341248, -4.11463999778592, 0.844031582886925, 7.78208960803073,
    11.5013829074586, 0.444482991372174, 2.60054581272466, -4.71288469356681,
    -0.343132143159804, -1.8798403907381, -0.31925032322771, 29.8395137569748,
    -19.713638420806, 1.21023454704822, -3.82133586653843, -4.24230669162403,
    -0.128656917758934, -2.00026616352991, 6.17406046755328, 1.21988831091319,
    0.755490574422338, 0.341153746411443, 8.78070330719093, 2.11077907947828,
    2.18271775066059, 0.623856270897387, 0.69746258657272, 0.306765097693532,
    2.55475109024386, -0.758870349117033, 0.0230517135428443, -10.4483044187463,
    -3.5425536379801, 0.50628903293276, 0.961044532294975, -2.4269356831673,
    -0.237630234314165, 28.4268931248783, 0.0796388029697632, 3.39006704666293,
    3.25845855709643, -0.482276238935394, -0.333295967717068, -0.251561675179636,
    11.9818989977435, 1.02120121164184, 3.83595253525676, -3.42827553387767,
    0.758245951293386, -0.572650437441927, 5.44181455783698, 0.195666536063851,
    -0.881636775462007, 0.637037329444995, 0.538262677735294, -3.63871523960419,
    -7.71004107609423, -0.229305695444564, -0.236313430436775, -1.0504337579587,
    0.169853833455182, 0.163739676205661, 0.849631675947778, 0.642716740413246,
    -0.84777150021684, 7.7193228935991, 1.52318785486357, 0.0948639492225293,
    0.461073658874404, -0.985275916004354, 0.584933835520111, 0.0876227706251002,
    -0.180712310040016, -0.84466888192973, 0.791051199519339,
    -0.0788890393807877, 2.56065499267076, -4.25565199655433, 0.905820269380403,
    0.381107537168204, 0.298752423017296, 0.757372425012899, -2.81908892546257,
    2.10289129235165, 2.15181714712737, 0.191185490396214, 61.1307866050661,
    5.61323651988839, 1.62336698760599, 25.0126955955884, 4.89427533544789,
    0.198301743557516, -0.926804721539349, 1.133842040068, 12.0824486879616,
    3.45221597386826, 0.640487320179048, -2.74368075373982, 0.577295692008625,
    -1.78206710207222, 2.5022615333175, 42.0633255270404, -1.53848460750582,
    24.9973838261109, 12.8816678499012, -2.82324502985904, 38.994311103535,
    2.02300868503141, 13.4151212631644, -2.49102773599107, -2.2922412046186,
    0.773021026947101, 0.891864337718481, 4.14965340394662, 0.127941799942574,
    26.6612202747941, 3.79007943734286, 82.5203173480609, 0.521486482945213,
    65.927083911124, 26.2222495778569, 1.43671637371689, -6.05869374608648,
    0.154614795278812, 14.9145956767882, -2.84998168439721, 3.79576005618016,
    -4.87858228794758, -2.3242386965776, 70.8000649505987, -0.0165902319188251,
    38.1861498961801, 42.5818268722957, 1.66543149866398, 4.19523710380519,
    -2.70888587615764, -2.60897826157812, 85.6839036318491, 9.51409185044254,
    2.07595242941192, 2.16914028993584, 425.282488345516, 2665.06275448215,
    983.038560966493, 1458.67332771568, 404.709088872717, 430.408747501351,
    146.56079353058, 183.125569821843, 233.279800191006, 134.580585125145,
    413.186901308111, 159.7416483475, 251.976931049503, 198.084433190346,
    357.597358409743, 1075.67503982347, 2.3444374746967, 2.29058228024483,
    2.03480625956375, 2.33712667685105, 3.35010888271603, 2.9007288258883,
    431.442383132788, 1190.94992978223, 418.57369363393, 558.866288339853,
    440.351570368775, 172.223318459985, 1343.04930401315, 571.692069897388,
    71.5262238106421, 44.4328526327163, 31.3754240863702, 18.3449028199435,
    194.760486764119, 2507.45453818833, 2715.77605824554, 218.241374434892,
    53.2229027089017, 27.0082886100378, 27.349542419444, 17.2289354255514,
    23.123001252027, 27.3152611148763, 27.2809557767091, 22.6640433474782,
    25.3097732134673, 35.9994765562133, 107.38995771048, 31.8636907211188,
    26.4551107803418, 196.671919429531, 197.237448621115, 18.0631255959236,
    27.177193258809, 40.8760525789647, 228.155127235726, 75.9641264577084,
    270.589618361396, 120.075948337205, 128.90935792457, 219.34183876427,
    107.964876862554, 46.8897452981378, 55.3345338921796, 39.0606798905224,
    37.5951737987281, 130.55818035623, 38.8976398275006, 77.0134404718054,
    33.7202040788422, 87.1624199438263, 257.706144840081, 299.136186572849,
    507.619557959078, 107.164352082557, 103.71537100687, 453.083762872283,
    365.318007037512, 2347.10841916393, 2127.89339954696, 1606.74241301685,
    314.379028205181, 169.470131293734, 296.389608036925, 120.745739747635,
    688.264946997008, 3016.64944380476, 3570.06249691295, 2989.5894034143,
    2892.54566145374, 2848.08615711965, 2868.18551172827, 2682.55554076977,
    2590.54279817676, 2577.28030598456, 2682.59062062179, 2968.17000914031,
    2775.7944333896, 2711.24321814979, 2813.37436349353, 2937.72100131844,
    3065.57416528862, 2713.41187156717, 2852.83112264489, 2936.92077309001,
    2802.92457204827, 2641.45699171615, 2751.39463474219, 2819.09871854182,
    4013.33170046984, 3102.88663336085, 2881.10573661047, 2743.52280434645,
    2767.49055656496, 3268.0691309758, 3317.43726120499, 3408.77613542638,
    3344.68547062379, 3249.42188912368, 2984.87543616333, 3096.31052304792,
    2820.30038469166, 2790.18009105722, 3019.82259564065, 2789.80683180017,
    2796.95728226352, 2655.16917653637, 3057.99288852393, 2862.31762501364,
    3662.041868057, 3507.62062063881, 3034.74683850852, 3405.79373852208,
    3048.72203391303, 2978.87467811178, 3088.68727275069, 3132.8427651132,
    3126.3518508079, 3252.29660363286, 3088.39243730773, 3009.0661041971,
    3215.6604385041, 2890.19467735806, 3074.96534773232, 2878.56122577563,
    2924.08604037177, 2934.66286812238, 2994.98722316008, 3025.24786218402,
    2896.50813109082, 3104.66999593572, 2979.96144691325, 3108.59514775849,
    3211.83438213981, 2919.00916604792, 2832.43146712377, 3354.79196850009,
    2868.02063652829, 2903.61398261496, 2841.48318102357, 2921.33456647919,
    2847.10170769894, 3824.86798809842, 3026.38311696209, 3174.03718674564,
    2933.88748913696, 2974.86577041583, 3024.06518561316, 2796.61745233139,
    2933.91870221052, 2870.66252320107, 2842.63917369635, 2865.2695615582,
    2958.41724008204, 2940.88916722871, 2995.14035142374, 3009.44095861494,
    2974.33404351043, 3005.93920501421, 2874.98216886936, 3191.22588375595,
    2951.46652485791, 3040.66338047994, 2967.5861580319, 3433.43654508035,
    5449.5265053144, 6007.74116414778, 5560.22310324003, 6023.81632874743,
    5687.98975280563, 5922.34770884077, 5992.15594524477, 5877.48005643479,
    5798.09878628578, 6194.20810440774, 5668.92071461122, 6564.46098726243,
    6145.56196475167, 5736.72147254813, 6730.36465122597, 6161.49389718292,
    6106.55685965108, 5827.65067993684, 5904.94077512375, 5798.16017137596,
    5911.16511118959, 5881.44766392918, 6281.52160352701, 5887.962208785,
    5796.71080252647, 6233.74191044089, 6404.03770809249, 5922.11945397977,
    5966.2755468527, 6041.22180032771, 6047.33934140136, 6151.3413167094,
    6739.55691876299, 6078.18642047924, 6050.05547319512, 6310.80631005949,
    6783.51996013186, 6657.54258862595, 6228.10480001357, 6549.22401573623,
    6531.91397268558, 6963.38352658672, 5666.1785104716, 6553.17812232908,
    7026.07015633034, 7212.44980180861, 6938.84523519744, 6709.37615308896,
    6535.44186277379, 6298.16502254648, 6650.67733667199, 6722.06641648871,
    6491.04640908057, 7008.48271196778, 9397.50723506215, 8217.28670302786,
    10471.0239182371, 0.0140692719979637, 0.0165774505746074, 31.7823672287145,
    62.9604028562611, 15.1664125413026, 26.2176569593046, 11.6093749867886,
    9.34483672055915, 4.43561288011778, 3.65704884593165, 2.90966929609479,
    1.91405052418481, 6.93196504220721, 3.06874984955933, 5.1647839141089,
    5.20946688251611, 7.11937629766632, 47.8446451070075, 0.023474345966044,
    0.0205663293304047, 0.0177581998239449, 0.0175839500174309,
    0.0793299419934655, 0.0642529255790744, 25.1296404255454, 30.0112074894203,
    6.3382153053485, 8.50333634593236, 9.93110892674407, 1.81330498512621,
    21.4401071418213, 13.6464333085191, 1.09522904433568, 0.72130240947977,
    0.559426746840054, 0.20559459931753, 8.07119017175418, 42.4966411061805,
    80.7461326667292, 5.16686076194696, 1.05278558147756, 0.293222647385074,
    0.404070260886554, 0.217744684034853, 0.23248468482683, 0.355670388338654,
    0.436233501004289, 0.198008471827508, 0.347305877870358, 0.693832311455581,
    2.26224052044736, 0.380902962958891, 0.358264617842252, 4.91273222676862,
    2.79990976060384, 0.215807182770692, 0.47787438206514, 0.887523284958117,
    6.06581463254804, 2.55586317369404, 4.11256886504178, 3.7930664482155,
    3.42128680823583, 3.29215850017993, 1.98056863091583, 0.504515052665318,
    1.01483861296915, 0.665105043431548, 0.459183812009783, 3.53307628712705,
    0.600058577518417, 1.93970239218434, 0.708504528232997, 3.69757250788691,
    6.62281963279432, 6.56378093865274, 14.3188773766095, 3.13161425520563,
    2.21170495436268, 9.09871256705425, 8.39104660870668, 84.8998740367256,
    30.3448314849108, 29.6327985268846, 6.59138553179582, 2.50441997576002,
    5.26310869792625, 1.5956748995866, 16.1910393320902, 22.843479727671,
    38.1521509414306, 20.3469452140349, 16.4168990720758, 18.0910158063208,
    18.1572204012465, 15.8259509777715, 17.6091088506813, 15.7971319111275,
    16.6880879650336, 36.6948750515511, 18.5114381083688, 17.6047370340549,
    17.7803571599497, 24.7595948331237, 28.7309741665407, 22.739827600408,
    23.2147360519811, 34.6928507574436, 20.1434241114332, 22.8746451251592,
    20.4645497628567, 19.3146276629425, 45.2006661067578, 22.2944941118367,
    19.8548163452791, 21.1551772770437, 17.8971553864412, 18.7965432642587,
    18.4918048867069, 18.5408883335957, 21.0867176894403, 29.59346063021,
    17.8537009606139, 19.7660101263783, 17.0432781665585, 17.8170071695916,
    19.0114183332166, 16.9621356321565, 14.1891632382502, 16.4354353709523,
    25.6195443890018, 16.5027354646262, 41.2531569913119, 37.617446004726,
    24.0042965549325, 31.8665296849699, 18.7964152896889, 21.7581073701846,
    21.1172489662034, 21.6232706685886, 21.8330240053743, 21.9025462287953,
    20.2138869146781, 20.1751944771102, 23.9924347547238, 20.6450484691328,
    27.3526061122805, 17.0679496130025, 19.2778062299836, 17.7796496343441,
    20.6439666942023, 21.2908146612719, 21.016951333208, 24.845636904401,
    19.9483009701992, 22.0140900174139, 21.8233603126833, 19.3524409535947,
    21.9534329953241, 35.7503186306607, 19.1952838008479, 18.3810956904515,
    19.6666091726204, 16.6317833799173, 16.5188297356357, 34.8066178001324,
    21.108655591925, 22.8488850983866, 21.9809714041449, 18.3838832392725,
    20.1105632575591, 17.747564435537, 17.6646987885885, 19.7287317904234,
    18.7449044421779, 18.4615048596345, 17.5043523232042, 15.6376841307986,
    19.450605226335, 19.73442835862, 21.1048050945396, 19.4348934949241,
    18.9134641238914, 19.5209441498156, 19.8731853326021, 19.536009113402,
    15.926733961238, 35.8631284065005, 55.9970670053903, 54.8900963575454,
    56.2814088936601, 59.5518624784346, 51.9015625211575, 57.2482145478397,
    54.7667490254746, 57.999913516774, 55.8336930603266, 52.4528628389068,
    54.024998010674, 67.6345360182024, 61.0099446737748, 50.9961762995641,
    61.954093114108, 56.5160037382618, 53.9437537128901, 51.962776732301,
    59.4962095973221, 51.1719152385972, 51.9640561375988, 53.455344964846,
    61.2788040479635, 53.8651913309016, 55.7966239334108, 52.606442167589,
    61.137030993731, 57.4820858335755, 59.6087871637874, 53.7154598845573,
    57.7153997710075, 57.1268591734244, 73.064784754887, 50.9585138759708,
    57.4858699024572, 60.5796366171701, 55.2360995194323, 63.1102560385421,
    57.2307219929518, 59.6345231601836, 59.3425695141662, 67.3282497973285,
    50.1625416754423, 59.3333511611789, 58.7073932847629, 50.6722669066123,
    59.0328168498566, 62.4136091104135, 52.1986438368176, 58.0015531360023,
    60.4986289746816, 58.5009413076356, 55.2176884072412, 66.7777562650089,
    77.8623493440223, 66.2746239641679, 82.4820183938641, 0.070729026815086,
    0.0963001482876482, 183.667683722993, 263.059292161335, 69.6618086537173,
    109.83685614149, 53.8372287354219, 48.1653858695095, 22.6313718693438,
    15.9379634078045, 13.8076416002474, 7.64195602681389, 30.7635875316098,
    18.8285528461362, 27.8854808805215, 31.9056069534251, 30.1479776593568,
    220.004747272928, 0.115091205344152, 0.090104401357782, 0.0778095428466797,
    0.0860429098970526, 0.352059496884784, 0.333485794161074, 115.840542385717,
    145.226727836446, 25.2159063662808, 33.491458322574, 40.0202003561212,
    10.1484866100891, 113.451543323037, 42.1571172853988, 5.0237575161657,
    4.86852890250687, 3.9314540376359, 0.843859124023288, 47.7979264373811,
    208.629352513027, 407.573821532521, 28.6341582139785, 4.50880830938839,
    1.34557746229254, 1.81323811676982, 1.23641467791307, 1.03127792908034,
    1.78494544127955, 2.36153300521655, 0.812396404392099, 1.51321047081691,
    3.40063593340756, 10.06059196214, 1.88577905844669, 1.72922495464854,
    18.1423967704437, 10.2953988265031, 1.12438473816713, 2.68334856580306,
    3.88309901567123, 25.3368551034608, 14.5055959740267, 19.9102119070044,
    20.6999612831075, 18.3719020532459, 14.211363121936, 9.5478282826769,
    2.29812283827159, 5.5911124909108, 3.00059163117087, 2.39627103268214,
    18.2249984413125, 3.56887951215307, 8.57392540382731, 3.9042438659397,
    21.7434278172237, 25.8582782626388, 27.7727656812938, 79.6114617672925,
    16.2226584601798, 10.7756072495855, 35.2994876794815, 37.178628709254,
    321.096509034326, 149.851321008285, 140.921088646662, 39.653558141417,
    12.8052091738065, 20.524587842587, 8.72815342144027, 76.7554947532622,
    108.48519934321, 186.099151165014, 103.045598784182, 88.9207019262156,
    90.3756273853443, 104.219719480106, 85.9437041310207, 82.8793260548144,
    90.0876682853601, 89.0285359970283, 175.090346210555, 86.9635767492356,
    88.30063252818, 90.6575089390534, 120.7670621762, 143.958406883025,
    107.056802069408, 99.6087395647432, 159.479673631948, 106.529310005859,
    102.574231415525, 93.5134804861962, 95.8574252629471, 194.675058560104,
    112.402736491307, 91.7871080595231, 95.6063432077198, 87.017572845427,
    101.041224515586, 108.680296672799, 122.181502737081, 118.536695554705,
    145.403867768185, 98.0295004938369, 103.085902640252, 104.037093236934,
    98.556048122494, 103.297709812703, 96.2476907036845, 88.6663845421293,
    86.667700742969, 123.406602561894, 97.7249529493807, 195.479459462024,
    182.400438440933, 106.390842417092, 154.925515761182, 99.8870132658013,
    126.279251378497, 99.1110855929892, 105.90241100353, 108.595271353368,
    123.600391186627, 99.0978594015137, 101.480807858051, 116.510336747994,
    104.805989369068, 131.047112860628, 89.4134934950644, 94.6444435664737,
    104.123380236309, 103.33047938675, 99.0055907599922, 102.174183508569,
    136.947131078562, 106.585863825735, 116.367249292971, 108.936717856644,
    92.3296205384942, 94.0932510838127, 170.525968279331, 94.0381261851433,
    104.268109746422, 97.8987232990514, 88.7846060482532, 92.7083938538146,
    158.189285646416, 125.178733017782, 111.396099389286, 111.866312722207,
    98.320324239371, 103.727782896833, 97.0225615428327, 87.9494375493342,
    98.9790194289132, 85.6656599598151, 97.5170224864035, 100.583448117768,
    99.2372247714617, 111.046218861364, 102.086928641283, 101.780714384496,
    107.086251177856, 98.6968956269254, 104.152771347386, 97.6633675119408,
    108.245079414401, 91.0733616082668, 184.85732245815, 237.904997328514,
    317.518534621323, 237.876543484401, 276.225873350593, 205.365016070203,
    247.581863801763, 269.497175497448, 282.242619353189, 260.745104807978,
    232.604288695607, 237.057901557766, 325.106260948821, 273.086067612923,
    252.725469593053, 273.391742540441, 263.465251491964, 251.28458266616,
    245.027030846171, 253.942544741242, 247.226515841904, 246.671083004086,
    239.3569827142, 262.697626361771, 236.962601589472, 225.418102190448,
    247.047570258848, 239.573909905756, 268.449969008803, 287.617289066228,
    227.054299698048, 236.506949364154, 257.022449648307, 382.954017415674,
    243.198743977303, 239.16374186172, 247.994023208933, 253.404525978926,
    276.628183167415, 230.958604175822, 235.487117454409, 259.839238157584,
    355.858859375677, 202.845551569225, 256.932265138715, 274.980085763363,
    247.889250928221, 285.245158284002, 280.114280556324, 266.618251711817,
    291.243651491869, 248.085239465447, 234.385178969872, 273.083371346997,
    275.964898607744, 338.026625495083, 341.145348752792, 356.39008273605,
    0.0111618622691932, 0.0176159420865874, 0.0176159420865874, 4.77123807850829,
    1.37970133279234, 1.96177440274405, 0.902730937667422, 0.969137508255345,
    0.468496490321149, 0.46028844176517, 0.779075657243679, 0.344281551937628,
    0.794489541860904, 0.285847067385881, 0.347638622121872, 0.483478559137166,
    0.778630608648657, 0.0176159420865874, 0.0111618622691932,
    0.0170873430500024, 0.015167922635657, 0.015167922635657, 0.0192735179083406,
    0.0170873430500024, 0.015167922635657, 1.07791143513693, 1.12458083867536,
    0.69543607914654, 0.482737332310643, 0.735077641719681, 3.815683331097,
    1.01173636690555, 0.244861189206639, 0.136994768504204, 0.0543030182093323,
    0.0836785794033855, 0.120781471621082, 5.74668072149755, 9.87521543545077,
    0.310487018291953, 0.100224428501467, 0.13820294435546, 0.026546564034222,
    0.015167922635657, 0.0922504440998582, 0.0578734061617651, 0.134911792654056,
    0.0678611761903362, 0.0420011091573045, 0.116221294527276, 0.251627243200058,
    0.109370457305903, 0.132129689366186, 0.0323595997043804, 0.218118048846752,
    0.0750086356311953, 0.0731932256469605, 0.10931803078345, 0.129582991449823,
    0.185397682399268, 0.901816042629452, 0.239913542522578, 0.142196377153014,
    0.36887332922921, 0.355688159177815, 0.183795273795009, 0.191202271535992,
    0.136802146341489, 0.10233377384197, 0.157304460746863, 0.114941682242679,
    0.105235450190731, 0.0716790286066608, 0.0087862497499484, 0.155269970600963,
    0.329117342526045, 0.413915090245805, 0.345185055949941, 0.330102882211229,
    0.399569496201098, 0.374146404463789, 0.62768661741852, 6.51473304340834,
    3.43064575022084, 0.779778286903276, 0.827836181942285, 0.768184056914801,
    0.302658459300356, 0.417763990085691, 19.3960426179089, 7.16921465239691,
    16.5105237267582, 18.8255507341192, 11.4178152202713, 19.6545685696086,
    22.2397013327665, 15.654652528526, 11.1250759051128, 23.9718760255547,
    10.8135739013821, 12.3156522367578, 17.5460220883701, 16.8036371634167,
    17.0759083506427, 14.7300892457821, 11.2686345117311, 12.8907589391191,
    13.1298510946008, 20.7810305767188, 10.3699987484981, 14.1513328722974,
    20.9512246996617, 21.5702344368644, 26.3337745517416, 23.9626813802805,
    7.76312991387212, 15.6684993519545, 19.0818445177699, 29.5433587912835,
    38.7788657518134, 29.6001264460249, 17.6869488359821, 28.9199303447177,
    17.1078630279312, 17.7570170721597, 21.8136306936345, 19.3173911701502,
    24.5316774744599, 23.9688314715559, 24.1861232359599, 17.8257356865846,
    27.6147651308716, 17.3562048740989, 26.0173021432172, 21.5425015082859,
    21.2000559668805, 21.758006749961, 17.2298748197572, 17.0069153186993,
    18.7790457529679, 28.8736268886103, 14.2562730643094, 20.9464040962439,
    16.0782788541495, 25.2023268461462, 16.6437823047442, 6.97161568181279,
    22.9796398992092, 20.8958644807413, 31.5867617870384, 14.7063637255257,
    14.2144106856554, 20.3906890544799, 23.2979257728298, 24.2527576061559,
    14.2860012386489, 12.7087418595287, 18.1978777956643, 16.3947779519285,
    17.8568399601377, 13.989181786293, 26.4578872514446, 16.4460683958068,
    18.703660995147, 18.5745948492296, 14.7174468863994, 19.875832967698,
    19.0607258905573, 17.134553811936, 17.4543421069913, 20.7492488132444,
    15.8417971679269, 26.128636483612, 24.9020615303017, 12.6022349372094,
    22.4280242581486, 23.593159301699, 34.4073502612392, 25.5083540149661,
    13.689201150752, 8.54237587514383, 24.2924622941496, 19.5751962888089,
    29.6248404873801, 27.8052743295075, 21.362135881168, 29.5514006941861,
    5.73070373913369, 27.5826595067271, 21.4039572358853, 11.1516570998728,
    24.6297016608789, 15.3496937135761, 25.8147984990087, 9.46557126517474,
    5.34873994057652, 31.3414425011795, 17.1962976477349, 16.7030990037076,
    37.0657872615075, 38.6250113065131, 21.7874517169151, 25.1594343729157,
    34.8621910427774, 29.6176775444156, 24.9398211997257, 5.21426440119382,
    21.3227370483483, 41.7113530153373, 9.38780151643863, 13.9366500744284,
    26.7847851504469, 9.8587342436645, 24.2531527446274, 17.3235979173201,
    33.7344803709299, 8.99080959270587, 28.3597175956683, 35.3247467393596,
    31.4097872595122, 8.89851989069157, 27.7246401589266, 10.4360866061152,
    10.4148833405425, 55.7884575978687, 12.29356410851, 6.15465057705031,
    24.4539768645127, 40.2024233189736, 45.6174421235178, 19.6095255213022,
    5.64996952285148, 16.5764821486622, 68.772142520343, 48.1058587973574,
    19.9845863908289, 56.6843761644099, 21.6164756225641, 21.7206236028976,
    23.0999199890822, 39.4916470387583, 4.69692878504838, 24.7581879122068,
    29.1665592049041, 70.6422766538908, 2.04544729146863, 2.12362840627665,
    381.789381533445, 2599.85758926462, 938.38557497397, 1439.12119001825,
    397.028695052893, 429.1345276559, 145.436098298974, 179.341852172797,
    227.715571741156, 131.202669452006, 399.827669494827, 155.614747718854,
    248.512849056472, 184.876983964979, 354.207372968382, 1074.30694717668,
    2.28875207454256, 2.2425540191095, 1.98402216645509, 2.27613873164915,
    3.29051247224626, 2.84914741289439, 387.641364451919, 1144.47347070554,
    411.184603778873, 557.075191648892, 439.306200801738, 168.975477340723,
    1304.11171788764, 553.030637911666, 69.8446273751714, 44.0254472218945,
    30.8266612313205, 17.7818552137251, 181.900264269517, 2485.77502111123,
    2699.0712180981, 214.732830248379, 52.6203651342827, 26.4700802452639,
    27.1027467198196, 17.0183630834116, 22.9072428013615, 26.7245334801598,
    26.1816368624113, 22.1224103804024, 24.4844600310535, 34.671791849133,
    105.748480724176, 31.3365772919353, 26.2512763604914, 193.779273391813,
    195.290988690399, 17.6948916466611, 26.747219650868, 40.0753190831744,
    227.585415516746, 68.4931403482612, 267.053765318234, 119.418200220258,
    120.593253958745, 213.89917806986, 107.266307577677, 46.3142541851819,
    54.3939164420536, 38.1907097069839, 37.046763316696, 129.299388767197,
    38.2375855373868, 76.3305490223698, 33.2436189461707, 86.302891640511,
    246.894450781451, 286.121206318071, 505.625944862419, 105.501427800572,
    102.895467129471, 448.618512498035, 359.240002309491, 2270.66856603252,
    2037.34651676596, 1587.34554285512, 298.485306705633, 164.338588912072,
    293.625327509079, 119.855932127374, 684.914895917272, 2960.11801686317,
    3498.65303150478, 2913.20342483514, 2826.72452359928, 2780.93280673631,
    2820.90062066785, 2638.39354409433, 2519.4925006874, 2543.1640459048,
    2630.35126461097, 2929.32309759165, 2735.98710668386, 2683.99331695963,
    2776.06363059757, 2860.79300549994, 2973.42887737047, 2647.99831908443,
    2776.02898098911, 2899.23300472439, 2730.32551384687, 2589.86945716674,
    2694.38695649368, 2751.11373471842, 3956.83234252247, 3048.51921450624,
    2849.29004962334, 2679.94785458677, 2733.69870268723, 3206.32954995844,
    3242.5962531448, 3337.07439607153, 3259.07316530547, 3188.31779621115,
    2915.7631834526, 3049.17199456877, 2767.23436252963, 2733.79385201763,
    2970.6526063392, 2733.7474176538, 2750.9796838407, 2597.4511183856,
    3000.42463999802, 2770.25873213043, 3573.63064475044, 3412.21750308176,
    2933.80096678675, 3319.21379493301, 2986.52841834106, 2910.97088731425,
    3005.94515916373, 3049.33163629658, 3058.85870044873, 3185.630267468,
    3005.62960586023, 2945.62078587481, 3134.93165163136, 2844.75405607522,
    3010.78416766181, 2833.95243394153, 2866.29157104748, 2878.57722044081,
    2940.05446832151, 2972.68778277824, 2848.55748874369, 3028.1033374308,
    2895.89106997235, 3077.48055274424, 3119.62998777812, 2876.38599495722,
    2786.5780616511, 3272.66732507914, 2784.89830521597, 2842.19812732374,
    2768.52196307069, 2869.2074537881, 2791.67720246186, 3717.85073300867,
    2984.44894999393, 3089.07657740689, 2856.42572156135, 2921.67552083911,
    2955.50729999663, 2750.00070757987, 2867.32013839514, 2807.13721077476,
    2776.17152712703, 2810.25955197015, 2914.70720132208, 2883.9610528403,
    2955.21975758144, 2947.91452776822, 2921.84710989191, 2939.93700541649,
    2836.17228469318, 3124.20949717208, 2866.42097024847, 2962.63348744327,
    2918.29461314047, 3381.79808375672, 5316.1044641904, 5878.13744349208,
    5479.66053744381, 5942.87118358882, 5598.88329073954, 5829.16277775156,
    5922.95846631311, 5823.98649597828, 5659.15947873271, 6086.29329120598,
    5472.77058638578, 6403.00336460982, 6080.14719471744, 5643.94078295452,
    6609.87352362962, 5984.07261365202, 5943.74714987836, 5670.50411477824,
    5745.19716539961, 5710.77509361552, 5770.43327013341, 5790.17642906898,
    6181.28219551626, 5841.45940858333, 5762.75267426247, 6135.10474982094,
    6346.58546461091, 5887.51125995469, 5879.42087483482, 5941.86653911728,
    5995.74291128569, 6093.03402289263, 6691.9351208196, 6002.06458213123,
    6006.42522733928, 6152.81974408482, 6622.67142962573, 6504.45577481254,
    6044.99416802702, 6456.16233770334, 6463.99872873035, 6910.77742027611,
    5637.19152950874, 6455.80144537152, 6864.49734957229, 7069.44713986063,
    6794.00922631744, 6522.10568257392, 6433.2754633599, 6183.39907004878,
    6494.69735529697, 6576.18664556264, 6377.65087878402, 6918.7474972295,
    9156.31454957305, 8025.29455399272, 10341.4525004054, 0.0, 0.0, 1.0, 5.0,
    0.0, 4.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 4.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 8.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
    0.0, 0.0, 3.0, 7.0, 5.0, 0.0, 0.0, 0.0, 0.0, 1.0, 15.0, 14.0, 16.0, 18.0,
    15.0, 14.0, 14.0, 17.0, 16.0, 16.0, 9.0, 14.0, 15.0, 15.0, 12.0, 10.0, 7.0,
    11.0, 8.0, 8.0, 7.0, 11.0, 14.0, 10.0, 16.0, 16.0, 15.0, 15.0, 15.0, 14.0,
    14.0, 13.0, 14.0, 15.0, 15.0, 17.0, 16.0, 15.0, 16.0, 14.0, 16.0, 14.0, 12.0,
    11.0, 15.0, 11.0, 13.0, 16.0, 13.0, 15.0, 16.0, 16.0, 15.0, 14.0, 15.0, 15.0,
    14.0, 13.0, 15.0, 16.0, 15.0, 17.0, 18.0, 16.0, 17.0, 13.0, 16.0, 17.0, 17.0,
    15.0, 13.0, 13.0, 15.0, 16.0, 15.0, 16.0, 13.0, 14.0, 15.0, 15.0, 17.0, 18.0,
    16.0, 16.0, 15.0, 18.0, 17.0, 17.0, 17.0, 14.0, 16.0, 14.0, 13.0, 16.0, 15.0,
    14.0, 16.0, 14.0, 12.0, 13.0, 13.0, 13.0, 11.0, 13.0, 12.0, 13.0, 13.0, 12.0,
    13.0, 13.0, 11.0, 12.0, 13.0, 12.0, 11.0, 14.0, 13.0, 14.0, 13.0, 13.0, 16.0,
    12.0, 12.0, 14.0, 14.0, 11.0, 14.0, 13.0, 12.0, 13.0, 12.0, 14.0, 12.0, 12.0,
    12.0, 14.0, 13.0, 14.0, 13.0, 14.0, 13.0, 14.0, 12.0, 12.0, 12.0, 12.0, 13.0,
    13.0, 12.0, 12.0, 14.0, 12.0, 12.0, 12.0, 9.0, 11.0, 0.0, 0.0, 0.0, 35.0,
    0.0, 28.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 26.0, 20.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 41.0, 20.0, 21.0, 0.0, 0.0, 0.0, 0.0, 0.0, 44.0, 45.0, 45.0, 45.0,
    45.0, 45.0, 46.0, 44.0, 45.0, 45.0, 42.0, 45.0, 44.0, 45.0, 38.0, 45.0, 39.0,
    42.0, 37.0, 38.0, 40.0, 46.0, 44.0, 44.0, 44.0, 45.0, 43.0, 43.0, 44.0, 42.0,
    45.0, 42.0, 47.0, 45.0, 45.0, 45.0, 44.0, 45.0, 46.0, 47.0, 41.0, 38.0, 38.0,
    40.0, 47.0, 35.0, 43.0, 46.0, 44.0, 45.0, 46.0, 45.0, 45.0, 44.0, 42.0, 44.0,
    43.0, 45.0, 46.0, 47.0, 46.0, 45.0, 47.0, 43.0, 46.0, 41.0, 46.0, 45.0, 46.0,
    40.0, 45.0, 40.0, 43.0, 44.0, 44.0, 44.0, 44.0, 42.0, 45.0, 45.0, 44.0, 46.0,
    47.0, 44.0, 44.0, 45.0, 46.0, 46.0, 44.0, 46.0, 46.0, 45.0, 43.0, 46.0, 47.0,
    43.0, 47.0, 44.0, 43.0, 45.0, 42.0, 45.0, 44.0, 44.0, 45.0, 46.0, 46.0, 44.0,
    45.0, 43.0, 43.0, 45.0, 44.0, 44.0, 44.0, 44.0, 44.0, 46.0, 41.0, 41.0, 47.0,
    46.0, 45.0, 45.0, 46.0, 45.0, 44.0, 44.0, 47.0, 46.0, 44.0, 42.0, 44.0, 43.0,
    43.0, 44.0, 44.0, 46.0, 45.0, 46.0, 44.0, 43.0, 45.0, 43.0, 43.0, 44.0, 44.0,
    44.0, 42.0, 47.0, 45.0, 45.0, 43.0, 44.0, 37.0, 45.0, 0.0, 0.0, 0.0, 7.0,
    0.0, 7.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    2.33333333333333, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.25, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3.25, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 13.6666666666667, 2.85714285714286, 4.2,
    0.0, 0.0, 0.0, 0.0, 0.0, 2.93333333333333, 3.21428571428571, 2.8125, 2.5,
    3.0, 3.21428571428571, 3.28571428571429, 2.58823529411765, 2.8125, 2.8125,
    4.66666666666667, 3.21428571428571, 2.93333333333333, 3.0, 3.16666666666667,
    4.5, 5.57142857142857, 3.81818181818182, 4.625, 4.75, 5.71428571428571,
    4.18181818181818, 3.14285714285714, 4.4, 2.75, 2.8125, 2.86666666666667,
    2.86666666666667, 2.93333333333333, 3.0, 3.21428571428571, 3.23076923076923,
    3.35714285714286, 3.0, 3.0, 2.64705882352941, 2.75, 3.0, 2.875,
    3.35714285714286, 2.5625, 2.71428571428571, 3.16666666666667,
    3.63636363636364, 3.13333333333333, 3.18181818181818, 3.30769230769231,
    2.875, 3.38461538461538, 3.0, 2.875, 2.8125, 3.0, 3.14285714285714, 2.8,
    2.93333333333333, 3.07142857142857, 3.46153846153846, 3.06666666666667,
    2.9375, 3.06666666666667, 2.64705882352941, 2.61111111111111, 2.6875,
    2.70588235294118, 3.15384615384615, 2.875, 2.64705882352941,
    2.70588235294118, 2.66666666666667, 3.46153846153846, 3.07692307692308,
    2.86666666666667, 2.75, 2.93333333333333, 2.75, 3.38461538461538, 3.0, 3.0,
    3.0, 2.58823529411765, 2.55555555555556, 2.9375, 2.75, 2.93333333333333, 2.5,
    2.70588235294118, 2.70588235294118, 2.58823529411765, 3.28571428571429,
    2.875, 3.21428571428571, 3.30769230769231, 2.875, 3.13333333333333,
    3.07142857142857, 2.9375, 3.14285714285714, 3.58333333333333,
    3.46153846153846, 3.23076923076923, 3.46153846153846, 4.0, 3.38461538461538,
    3.75, 3.53846153846154, 3.53846153846154, 3.66666666666667, 3.46153846153846,
    3.30769230769231, 3.90909090909091, 3.75, 3.38461538461538, 3.66666666666667,
    4.0, 3.14285714285714, 3.38461538461538, 3.28571428571429, 3.15384615384615,
    3.15384615384615, 2.9375, 3.83333333333333, 3.75, 3.21428571428571,
    3.28571428571429, 4.09090909090909, 3.14285714285714, 3.38461538461538,
    3.91666666666667, 3.53846153846154, 3.66666666666667, 3.0, 3.66666666666667,
    3.58333333333333, 3.58333333333333, 3.14285714285714, 3.38461538461538,
    3.28571428571429, 3.46153846153846, 3.28571428571429, 3.38461538461538,
    3.07142857142857, 3.75, 3.58333333333333, 3.58333333333333, 3.66666666666667,
    3.38461538461538, 3.38461538461538, 3.5, 3.91666666666667, 3.21428571428571,
    3.75, 3.58333333333333, 3.66666666666667, 4.11111111111111, 4.09090909090909,
    0.0, 0.0, 1.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 3.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 6.0, 5.0, 6.0, 8.0, 7.0, 7.0, 6.0,
    7.0, 6.0, 7.0, 7.0, 12.0, 7.0, 6.0, 12.0, 9.0, 7.0, 6.0, 7.0, 7.0, 7.0, 8.0,
    10.0, 7.0, 9.0, 11.0, 11.0, 6.0, 8.0, 8.0, 11.0, 9.0, 14.0, 7.0, 9.0, 11.0,
    12.0, 13.0, 14.0, 15.0, 13.0, 15.0, 9.0, 14.0, 17.0, 17.0, 15.0, 11.0, 11.0,
    10.0, 14.0, 15.0, 10.0, 14.0, 27.0, 20.0, 32.0, 0.0, 0.0, 1.0, 2.0, 0.0, 3.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 8.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    5.0, 4.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 8.0, 2.0, 0.0, 0.0, 3.0, 0.0,
    0.0, 0.0, 0.0, 6.0, 0.0, 0.0, 0.0, 5.0, 8.0, 4.0, 4.0, 8.0, 4.0, 3.0, 3.0,
    2.0, 14.0, 6.0, 0.0, 2.0, 0.0, 5.0, 5.0, 5.0, 6.0, 9.0, 3.0, 2.0, 1.0, 1.0,
    3.0, 1.0, 0.0, 0.0, 6.0, 3.0, 8.0, 5.0, 6.0, 10.0, 3.0, 4.0, 4.0, 5.0, 7.0,
    3.0, 4.0, 4.0, 9.0, 2.0, 5.0, 0.0, 2.0, 2.0, 3.0, 2.0, 3.0, 4.0, 2.0, 5.0,
    4.0, 2.0, 3.0, 9.0, 3.0, 1.0, 2.0, 0.0, 1.0, 14.0, 3.0, 6.0, 4.0, 2.0, 4.0,
    1.0, 0.0, 2.0, 0.0, 1.0, 2.0, 1.0, 4.0, 3.0, 5.0, 2.0, 2.0, 5.0, 2.0, 4.0,
    0.0, 10.0, 26.0, 31.0, 25.0, 25.0, 25.0, 24.0, 29.0, 24.0, 27.0, 27.0, 24.0,
    23.0, 26.0, 25.0, 23.0, 24.0, 27.0, 27.0, 26.0, 24.0, 24.0, 21.0, 23.0, 22.0,
    22.0, 26.0, 23.0, 24.0, 22.0, 25.0, 21.0, 25.0, 22.0, 28.0, 24.0, 23.0, 24.0,
    21.0, 20.0, 19.0, 21.0, 21.0, 22.0, 22.0, 21.0, 24.0, 20.0, 27.0, 24.0, 23.0,
    21.0, 22.0, 26.0, 25.0, 19.0, 25.0, 15.0, 50.0, 50.0, 48.0, 43.0, 50.0, 47.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 47.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 48.0, 49.0, 50.0, 50.0, 50.0, 50.0, 49.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 41.0, 44.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 40.0, 46.0, 46.0,
    50.0, 50.0, 50.0, 50.0, 50.0, 45.0, 40.0, 48.0, 50.0, 50.0, 47.0, 50.0, 50.0,
    50.0, 50.0, 43.0, 50.0, 50.0, 50.0, 45.0, 42.0, 46.0, 46.0, 42.0, 46.0, 47.0,
    47.0, 48.0, 34.0, 44.0, 50.0, 48.0, 50.0, 45.0, 45.0, 45.0, 44.0, 41.0, 47.0,
    48.0, 49.0, 49.0, 47.0, 49.0, 50.0, 50.0, 44.0, 47.0, 40.0, 42.0, 44.0, 40.0,
    47.0, 46.0, 46.0, 45.0, 43.0, 47.0, 46.0, 46.0, 41.0, 48.0, 45.0, 50.0, 48.0,
    48.0, 47.0, 48.0, 47.0, 46.0, 48.0, 45.0, 46.0, 48.0, 47.0, 41.0, 47.0, 49.0,
    48.0, 50.0, 49.0, 36.0, 47.0, 44.0, 46.0, 48.0, 46.0, 49.0, 50.0, 48.0, 50.0,
    49.0, 48.0, 49.0, 46.0, 47.0, 45.0, 48.0, 48.0, 45.0, 48.0, 46.0, 50.0, 39.0,
    18.0, 14.0, 19.0, 17.0, 18.0, 19.0, 15.0, 19.0, 17.0, 16.0, 19.0, 15.0, 17.0,
    19.0, 15.0, 17.0, 16.0, 17.0, 17.0, 19.0, 19.0, 21.0, 17.0, 21.0, 19.0, 13.0,
    16.0, 20.0, 20.0, 17.0, 18.0, 16.0, 14.0, 15.0, 17.0, 16.0, 14.0, 16.0, 16.0,
    16.0, 16.0, 14.0, 19.0, 14.0, 12.0, 9.0, 15.0, 12.0, 15.0, 17.0, 15.0, 13.0,
    14.0, 11.0, 4.0, 5.0, 3.0 };

  double dist[245];
  int c_k;
  int idx[245];
  double s;
  signed char b_Y[3];
  emxArray_real_T *edges;
  int exitg4;
  int exitg3;
  emxArray_real_T *nn;
  unsigned char outsize_idx_0;
  boolean_T guard1 = false;
  int exitg2;
  int high_i;
  int mid_i;
  boolean_T exitg1;
  emxInit_real_T1(&Uc, 1);
  low_ip1 = Uc->size[0];
  Uc->size[0] = 245;
  emxEnsureCapacity((emxArray__common *)Uc, low_ip1, (int)sizeof(double));
  for (k = 0; k < 245; k++) {
    Uc->data[k] = a[k];
  }

  nb = -1;
  k = 0;
  while (k + 1 <= 245) {
    nbins = (int)Uc->data[k];
    do {
      exitg5 = 0;
      k++;
      if (k + 1 > 245) {
        exitg5 = 1;
      } else {
        eok = (std::abs((double)nbins - Uc->data[k]) < eps((double)nbins / 2.0));
        if (!eok) {
          exitg5 = 1;
        }
      }
    } while (exitg5 == 0);

    nb++;
    Uc->data[nb] = nbins;
  }

  emxInit_real_T(&b_testSamples, 2);
  low_ip1 = Uc->size[0];
  if (1 > nb + 1) {
    i7 = -1;
  } else {
    i7 = nb;
  }

  Uc->size[0] = i7 + 1;
  emxEnsureCapacity((emxArray__common *)Uc, low_ip1, (int)sizeof(double));
  nbins = testSamples->size[1];
  low_ip1 = b_testSamples->size[0] * b_testSamples->size[1];
  b_testSamples->size[0] = 1;
  b_testSamples->size[1] = nbins;
  emxEnsureCapacity((emxArray__common *)b_testSamples, low_ip1, (int)sizeof
                    (double));
  for (low_ip1 = 0; low_ip1 < nbins; low_ip1++) {
    b_testSamples->data[b_testSamples->size[0] * low_ip1] = testSamples->
      data[testSamples->size[0] * low_ip1];
  }

  emxInit_real_T(&r4, 2);
  low_ip1 = r4->size[0] * r4->size[1];
  r4->size[0] = 245;
  r4->size[1] = b_testSamples->size[1];
  emxEnsureCapacity((emxArray__common *)r4, low_ip1, (int)sizeof(double));
  for (low_ip1 = 0; low_ip1 < 245; low_ip1++) {
    nbins = b_testSamples->size[1];
    for (nb = 0; nb < nbins; nb++) {
      r4->data[low_ip1 + r4->size[0] * nb] = b_testSamples->data
        [b_testSamples->size[0] * nb];
    }
  }

  emxFree_real_T(&b_testSamples);
  for (low_ip1 = 0; low_ip1 < 104; low_ip1++) {
    for (nb = 0; nb < 245; nb++) {
      x[nb + 245 * low_ip1] = tX[nb + 245 * low_ip1] - r4->data[nb + 245 *
        low_ip1];
    }
  }

  emxFree_real_T(&r4);

  for (b_k = 1; b_k < 25481; b_k++) {
    c_k = b_k;
    y[c_k - 1] = x[c_k - 1] * x[c_k - 1];
  }

  for (nbins = 0; nbins < 245; nbins++) {
    s = y[nbins];
    for (k = 0; k < 103; k++) {
      s += y[nbins + (k + 1) * 245];
    }

    dist[nbins] = s;
  }

  c_sort(dist, idx);
  for (nbins = 0; nbins < 3; nbins++) {
    b_Y[nbins] = a[idx[nbins] - 1];
  }

  emxInit_real_T(&edges, 2);
  nbins = Uc->size[0];
  low_ip1 = edges->size[0] * edges->size[1];
  edges->size[0] = 1;
  edges->size[1] = (unsigned char)(nbins + 1);
  emxEnsureCapacity((emxArray__common *)edges, low_ip1, (int)sizeof(double));
  k = 0;
  do {
    exitg4 = 0;
    nbins = Uc->size[0];
    if (k <= nbins - 2) {
      edges->data[1 + k] = Uc->data[k] + (Uc->data[1 + k] - Uc->data[k]) / 2.0;
      k++;
    } else {
      exitg4 = 1;
    }
  } while (exitg4 == 0);

  edges->data[0] = rtMinusInf;
  edges->data[edges->size[1] - 1] = rtInf;
  k = 1;
  do {
    exitg3 = 0;
    nbins = Uc->size[0];
    if (k - 1 <= nbins - 2) {
      edges->data[k] += eps(edges->data[k]);
      k++;
    } else {
      exitg3 = 1;
    }
  } while (exitg3 == 0);

  emxInit_real_T1(&nn, 1);
  outsize_idx_0 = (unsigned char)edges->size[1];
  low_ip1 = nn->size[0];
  nn->size[0] = outsize_idx_0;
  emxEnsureCapacity((emxArray__common *)nn, low_ip1, (int)sizeof(double));
  nbins = outsize_idx_0;
  for (low_ip1 = 0; low_ip1 < nbins; low_ip1++) {
    nn->data[low_ip1] = 0.0;
  }

  nbins = edges->size[1];
  guard1 = false;
  if (nbins > 1) {
    nb = 1;
    do {
      exitg2 = 0;
      if (nb + 1 <= nbins) {
        if (!(edges->data[nb] >= edges->data[nb - 1])) {
          eok = false;
          exitg2 = 1;
        } else {
          nb++;
        }
      } else {
        guard1 = true;
        exitg2 = 1;
      }
    } while (exitg2 == 0);
  } else {
    guard1 = true;
  }

  if (guard1) {
    eok = true;
  }

  if (!eok) {
    low_ip1 = nn->size[0];
    nn->size[0] = outsize_idx_0;
    emxEnsureCapacity((emxArray__common *)nn, low_ip1, (int)sizeof(double));
    nbins = outsize_idx_0;
    for (low_ip1 = 0; low_ip1 < nbins; low_ip1++) {
      nn->data[low_ip1] = rtNaN;
    }
  } else {
    nbins = 0;
    for (k = 0; k < 3; k++) {
      nb = 0;
      if ((b_Y[nbins] >= edges->data[0]) && (b_Y[nbins] < edges->data
           [edges->size[1] - 1])) {
        nb = 1;
        low_ip1 = 2;
        high_i = edges->size[1];
        while (high_i > low_ip1) {
          mid_i = (nb >> 1) + (high_i >> 1);
          if (((nb & 1) == 1) && ((high_i & 1) == 1)) {
            mid_i++;
          }

          if (b_Y[nbins] >= edges->data[mid_i - 1]) {
            nb = mid_i;
            low_ip1 = mid_i + 1;
          } else {
            high_i = mid_i;
          }
        }
      }

      if (b_Y[nbins] == edges->data[edges->size[1] - 1]) {
        nb = edges->size[1];
      }

      if (nb > 0) {
        nn->data[nb - 1]++;
      }

      nbins++;
    }
  }

  low_ip1 = edges->size[0] * edges->size[1];
  edges->size[0] = 1;
  edges->size[1] = nn->size[0] - 1;
  emxEnsureCapacity((emxArray__common *)edges, low_ip1, (int)sizeof(double));
  for (k = 0; k <= nn->size[0] - 2; k++) {
    edges->data[k] = nn->data[k];
  }

  if (nn->size[0] - 1 > 0) {
    edges->data[edges->size[1] - 1] += nn->data[nn->size[0] - 1];
  }

  emxFree_real_T(&nn);
  nbins = 1;
  nb = edges->size[1];
  s = edges->data[0];
  low_ip1 = 0;
  if (edges->size[1] > 1) {
    if (rtIsNaN(edges->data[0])) {
      high_i = 2;
      exitg1 = false;
      while ((!exitg1) && (high_i <= nb)) {
        nbins = high_i;
        if (!rtIsNaN(edges->data[high_i - 1])) {
          s = edges->data[high_i - 1];
          low_ip1 = high_i - 1;
          exitg1 = true;
        } else {
          high_i++;
        }
      }
    }

    if (nbins < edges->size[1]) {
      while (nbins + 1 <= nb) {
        if (edges->data[nbins] > s) {
          s = edges->data[nbins];
          low_ip1 = nbins;
        }

        nbins++;
      }
    }
  }

  emxFree_real_T(&edges);
  Y = Uc->data[low_ip1];
  emxFree_real_T(&Uc);
  return Y;
}

//
// Arguments    : const emxArray_real_T *x
// Return Type  : double
//
static double mean(const emxArray_real_T *x)
{
  double y;
  int k;
  if (x->size[0] == 0) {
    y = 0.0;
  } else {
    y = x->data[0];
    for (k = 2; k <= x->size[0]; k++) {
      y += x->data[k - 1];
    }
  }

  y /= (double)x->size[0];
  return y;
}

//
// Arguments    : emxArray_int32_T *idx
//                emxArray_real_T *x
//                int offset
//                int np
//                int nq
//                emxArray_int32_T *iwork
//                emxArray_real_T *xwork
// Return Type  : void
//
static void merge(emxArray_int32_T *idx, emxArray_real_T *x, int offset, int np,
                  int nq, emxArray_int32_T *iwork, emxArray_real_T *xwork)
{
  int n;
  int qend;
  int p;
  int iout;
  int exitg1;
  if ((np == 0) || (nq == 0)) {
  } else {
    n = np + nq;
    for (qend = 0; qend + 1 <= n; qend++) {
      iwork->data[qend] = idx->data[offset + qend];
      xwork->data[qend] = x->data[offset + qend];
    }

    p = 0;
    n = np;
    qend = np + nq;
    iout = offset - 1;
    do {
      exitg1 = 0;
      iout++;
      if (xwork->data[p] >= xwork->data[n]) {
        idx->data[iout] = iwork->data[p];
        x->data[iout] = xwork->data[p];
        if (p + 1 < np) {
          p++;
        } else {
          exitg1 = 1;
        }
      } else {
        idx->data[iout] = iwork->data[n];
        x->data[iout] = xwork->data[n];
        if (n + 1 < qend) {
          n++;
        } else {
          n = iout - p;
          while (p + 1 <= np) {
            idx->data[(n + p) + 1] = iwork->data[p];
            x->data[(n + p) + 1] = xwork->data[p];
            p++;
          }

          exitg1 = 1;
        }
      }
    } while (exitg1 == 0);
  }
}

//
// Arguments    : emxArray_int32_T *idx
//                emxArray_real_T *x
//                int offset
//                int n
//                int preSortLevel
//                emxArray_int32_T *iwork
//                emxArray_real_T *xwork
// Return Type  : void
//
static void merge_block(emxArray_int32_T *idx, emxArray_real_T *x, int offset,
  int n, int preSortLevel, emxArray_int32_T *iwork, emxArray_real_T *xwork)
{
  int nPairs;
  int bLen;
  int tailOffset;
  int nTail;
  nPairs = n >> preSortLevel;
  bLen = 1 << preSortLevel;
  while (nPairs > 1) {
    if ((nPairs & 1) != 0) {
      nPairs--;
      tailOffset = bLen * nPairs;
      nTail = n - tailOffset;
      if (nTail > bLen) {
        merge(idx, x, offset + tailOffset, bLen, nTail - bLen, iwork, xwork);
      }
    }

    tailOffset = bLen << 1;
    nPairs >>= 1;
    for (nTail = 1; nTail <= nPairs; nTail++) {
      merge(idx, x, offset + (nTail - 1) * tailOffset, bLen, bLen, iwork, xwork);
    }

    bLen = tailOffset;
  }

  if (n > bLen) {
    merge(idx, x, offset, bLen, n - bLen, iwork, xwork);
  }
}

//
// Arguments    : const emxArray_real_T *A
//                const emxArray_real_T *B
//                emxArray_real_T *y
// Return Type  : void
//
static void mrdivide(const emxArray_real_T *A, const emxArray_real_T *B,
                     emxArray_real_T *y)
{
  emxArray_real_T *Y;
  emxArray_real_T *b_B;
  emxArray_real_T *b_A;
  emxArray_real_T *tau;
  emxArray_int32_T *jpvt;
  emxArray_real_T *c_B;
  unsigned int unnamed_idx_0;
  int i2;
  double tol;
  unsigned int unnamed_idx_1;
  double temp;
  int minmn;
  int k;
  int rankR;
  int maxmn;
  emxInit_real_T(&Y, 2);
  emxInit_real_T(&b_B, 2);
  emxInit_real_T(&b_A, 2);
  emxInit_real_T1(&tau, 1);
  emxInit_int32_T1(&jpvt, 2);
  emxInit_real_T(&c_B, 2);
  if ((A->size[0] == 0) || (B->size[0] == 0)) {
    unnamed_idx_0 = (unsigned int)A->size[0];
    unnamed_idx_1 = (unsigned int)B->size[0];
    i2 = y->size[0] * y->size[1];
    y->size[0] = (int)unnamed_idx_0;
    y->size[1] = (int)unnamed_idx_1;
    emxEnsureCapacity((emxArray__common *)y, i2, (int)sizeof(double));
    k = (int)unnamed_idx_0 * (int)unnamed_idx_1;
    for (i2 = 0; i2 < k; i2++) {
      y->data[i2] = 0.0;
    }
  } else if (B->size[0] == 1) {
    tol = A->data[0];
    temp = 1.0 / B->data[0];
    minmn = 1;
    while (minmn <= A->size[0]) {
      tol *= temp;
      minmn = 2;
    }

    i2 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)y, i2, (int)sizeof(double));
    y->data[0] = tol;
  } else {
    i2 = b_B->size[0] * b_B->size[1];
    b_B->size[0] = 1;
    b_B->size[1] = A->size[0];
    emxEnsureCapacity((emxArray__common *)b_B, i2, (int)sizeof(double));
    k = A->size[0];
    for (i2 = 0; i2 < k; i2++) {
      b_B->data[b_B->size[0] * i2] = A->data[i2];
    }

    i2 = b_A->size[0] * b_A->size[1];
    b_A->size[0] = 1;
    b_A->size[1] = B->size[0];
    emxEnsureCapacity((emxArray__common *)b_A, i2, (int)sizeof(double));
    k = B->size[0];
    for (i2 = 0; i2 < k; i2++) {
      b_A->data[b_A->size[0] * i2] = B->data[i2];
    }

    xgeqp3(b_A, tau, jpvt);
    rankR = 0;
    if (1 < b_A->size[1]) {
      minmn = 1;
      maxmn = b_A->size[1];
    } else {
      minmn = b_A->size[1];
      maxmn = 1;
    }

    if (minmn > 0) {
      tol = (double)maxmn * std::abs(b_A->data[0]) * 2.2204460492503131E-16;
      while ((rankR < 1) && (std::abs(b_A->data[0]) >= tol)) {
        rankR = 1;
      }
    }

    minmn = b_A->size[1];
    maxmn = b_B->size[1];
    i2 = Y->size[0] * Y->size[1];
    Y->size[0] = minmn;
    Y->size[1] = maxmn;
    emxEnsureCapacity((emxArray__common *)Y, i2, (int)sizeof(double));
    k = minmn * maxmn;
    for (i2 = 0; i2 < k; i2++) {
      Y->data[i2] = 0.0;
    }

    i2 = c_B->size[0] * c_B->size[1];
    c_B->size[0] = 1;
    c_B->size[1] = b_B->size[1];
    emxEnsureCapacity((emxArray__common *)c_B, i2, (int)sizeof(double));
    k = b_B->size[0] * b_B->size[1];
    for (i2 = 0; i2 < k; i2++) {
      c_B->data[i2] = b_B->data[i2];
    }

    minmn = b_A->size[1];
    if (1 <= minmn) {
      minmn = 1;
    }

    maxmn = 1;
    while (maxmn <= minmn) {
      if (tau->data[0] != 0.0) {
        for (k = 0; k + 1 <= b_B->size[1]; k++) {
          tol = c_B->data[c_B->size[0] * k];
          tol *= tau->data[0];
          if (tol != 0.0) {
            c_B->data[c_B->size[0] * k] -= tol;
          }
        }
      }

      maxmn = 2;
    }

    for (k = 0; k + 1 <= b_B->size[1]; k++) {
      minmn = 1;
      while (minmn <= rankR) {
        Y->data[(jpvt->data[0] + Y->size[0] * k) - 1] = c_B->data[c_B->size[0] *
          k];
        minmn = 2;
      }

      maxmn = rankR;
      while (maxmn > 0) {
        Y->data[(jpvt->data[0] + Y->size[0] * k) - 1] /= b_A->data[0];
        maxmn = 0;
      }
    }

    i2 = y->size[0] * y->size[1];
    y->size[0] = Y->size[1];
    y->size[1] = Y->size[0];
    emxEnsureCapacity((emxArray__common *)y, i2, (int)sizeof(double));
    k = Y->size[0];
    for (i2 = 0; i2 < k; i2++) {
      minmn = Y->size[1];
      for (maxmn = 0; maxmn < minmn; maxmn++) {
        y->data[maxmn + y->size[0] * i2] = Y->data[i2 + Y->size[0] * maxmn];
      }
    }
  }

  emxFree_real_T(&c_B);
  emxFree_int32_T(&jpvt);
  emxFree_real_T(&tau);
  emxFree_real_T(&b_A);
  emxFree_real_T(&b_B);
  emxFree_real_T(&Y);
}

//
// Arguments    : const double Y[50]
//                emxArray_real_T *iPk
//                double Ph
// Return Type  : void
//
static void removePeaksBelowMinPeakHeight(const double Y[50], emxArray_real_T
  *iPk, double Ph)
{
  int end;
  int trueCount;
  int i;
  int partialTrueCount;
  if (!(iPk->size[0] == 0)) {
    end = iPk->size[0] - 1;
    trueCount = 0;
    for (i = 0; i <= end; i++) {
      if (Y[(int)iPk->data[i] - 1] > Ph) {
        trueCount++;
      }
    }

    partialTrueCount = 0;
    for (i = 0; i <= end; i++) {
      if (Y[(int)iPk->data[i] - 1] > Ph) {
        iPk->data[partialTrueCount] = iPk->data[i];
        partialTrueCount++;
      }
    }

    end = iPk->size[0];
    iPk->size[0] = trueCount;
    emxEnsureCapacity((emxArray__common *)iPk, end, (int)sizeof(double));
  }
}

//
// Arguments    : const double Y[50]
//                emxArray_real_T *iPk
//                double Th
// Return Type  : void
//
static void removePeaksBelowThreshold(const double Y[50], emxArray_real_T *iPk,
  double Th)
{
  int c;
  emxArray_real_T *base;
  int k;
  int trueCount;
  double extremum;
  int partialTrueCount;
  c = iPk->size[0];
  emxInit_real_T1(&base, 1);
  k = base->size[0];
  base->size[0] = c;
  emxEnsureCapacity((emxArray__common *)base, k, (int)sizeof(double));
  for (k = 0; k + 1 <= c; k++) {
    if ((Y[(int)(iPk->data[k] - 1.0) - 1] >= Y[(int)(iPk->data[k] + 1.0) - 1]) ||
        rtIsNaN(Y[(int)(iPk->data[k] + 1.0) - 1])) {
      extremum = Y[(int)(iPk->data[k] - 1.0) - 1];
    } else {
      extremum = Y[(int)(iPk->data[k] + 1.0) - 1];
    }

    base->data[k] = extremum;
  }

  k = iPk->size[0] - 1;
  trueCount = 0;
  for (c = 0; c <= k; c++) {
    if (Y[(int)iPk->data[c] - 1] - base->data[c] >= Th) {
      trueCount++;
    }
  }

  partialTrueCount = 0;
  for (c = 0; c <= k; c++) {
    if (Y[(int)iPk->data[c] - 1] - base->data[c] >= Th) {
      iPk->data[partialTrueCount] = iPk->data[c];
      partialTrueCount++;
    }
  }

  emxFree_real_T(&base);
  k = iPk->size[0];
  iPk->size[0] = trueCount;
  emxEnsureCapacity((emxArray__common *)iPk, k, (int)sizeof(double));
}

//
// Get Resultant. Size of input vector must be greater than 3.
// Arguments    : const double X[150]
//                double Y[50]
// Return Type  : void
//
static void resultant(const double X[150], double Y[50])
{
  double X1[150];
  int k;
  for (k = 0; k < 150; k++) {
    X1[k] = X[k] * X[k];
  }

  for (k = 0; k < 50; k++) {
    Y[k] = std::sqrt((X1[k] + X1[50 + k]) + X1[100 + k]);
  }
}

//
// Arguments    : double u0
//                double u1
// Return Type  : double
//
static double rt_hypotd_snf(double u0, double u1)
{
  double y;
  double a;
  double b;
  a = std::abs(u0);
  b = std::abs(u1);
  if (a < b) {
    a /= b;
    y = b * std::sqrt(a * a + 1.0);
  } else if (a > b) {
    b /= a;
    y = a * std::sqrt(b * b + 1.0);
  } else if (rtIsNaN(b)) {
    y = b;
  } else {
    y = a * 1.4142135623730951;
  }

  return y;
}

//
// Arguments    : int *k
//                const emxArray_real_T *x
// Return Type  : double
//
static double skip_to_last_equal_value(int *k, const emxArray_real_T *x)
{
  double xk;
  boolean_T exitg1;
  boolean_T p;
  xk = x->data[*k - 1];
  exitg1 = false;
  while ((!exitg1) && (*k < x->size[0])) {
    if ((std::abs(xk - x->data[*k]) < eps(xk / 2.0)) || (rtIsInf(x->data[*k]) &&
         rtIsInf(xk) && ((x->data[*k] > 0.0) == (xk > 0.0)))) {
      p = true;
    } else {
      p = false;
    }

    if (p) {
      (*k)++;
    } else {
      exitg1 = true;
    }
  }

  return xk;
}

//
// Arguments    : emxArray_real_T *x
//                emxArray_int32_T *idx
// Return Type  : void
//
static void sort(emxArray_real_T *x, emxArray_int32_T *idx)
{
  int dim;
  dim = 2;
  if (x->size[0] != 1) {
    dim = 1;
  }

  b_sort(x, dim, idx);
}

//
// Arguments    : emxArray_real_T *x
//                emxArray_int32_T *idx
// Return Type  : void
//
static void sortIdx(emxArray_real_T *x, emxArray_int32_T *idx)
{
  emxArray_real_T *b_x;
  unsigned int unnamed_idx_0;
  int ib;
  int m;
  int n;
  double x4[4];
  int idx4[4];
  emxArray_int32_T *iwork;
  emxArray_real_T *xwork;
  int nNaNs;
  int k;
  int wOffset;
  signed char perm[4];
  int nNonNaN;
  int p;
  int i4;
  int nBlocks;
  int b_iwork[256];
  double b_xwork[256];
  int b;
  int bLen;
  int bLen2;
  int nPairs;
  int exitg1;
  emxInit_real_T1(&b_x, 1);
  unnamed_idx_0 = (unsigned int)x->size[0];
  ib = b_x->size[0];
  b_x->size[0] = x->size[0];
  emxEnsureCapacity((emxArray__common *)b_x, ib, (int)sizeof(double));
  m = x->size[0];
  for (ib = 0; ib < m; ib++) {
    b_x->data[ib] = x->data[ib];
  }

  ib = idx->size[0];
  idx->size[0] = (int)unnamed_idx_0;
  emxEnsureCapacity((emxArray__common *)idx, ib, (int)sizeof(int));
  m = (int)unnamed_idx_0;
  for (ib = 0; ib < m; ib++) {
    idx->data[ib] = 0;
  }

  n = x->size[0];
  for (m = 0; m < 4; m++) {
    x4[m] = 0.0;
    idx4[m] = 0;
  }

  emxInit_int32_T(&iwork, 1);
  ib = iwork->size[0];
  iwork->size[0] = (int)unnamed_idx_0;
  emxEnsureCapacity((emxArray__common *)iwork, ib, (int)sizeof(int));
  m = iwork->size[0];
  ib = iwork->size[0];
  iwork->size[0] = m;
  emxEnsureCapacity((emxArray__common *)iwork, ib, (int)sizeof(int));
  for (ib = 0; ib < m; ib++) {
    iwork->data[ib] = 0;
  }

  emxInit_real_T1(&xwork, 1);
  unnamed_idx_0 = (unsigned int)x->size[0];
  ib = xwork->size[0];
  xwork->size[0] = (int)unnamed_idx_0;
  emxEnsureCapacity((emxArray__common *)xwork, ib, (int)sizeof(double));
  m = xwork->size[0];
  ib = xwork->size[0];
  xwork->size[0] = m;
  emxEnsureCapacity((emxArray__common *)xwork, ib, (int)sizeof(double));
  for (ib = 0; ib < m; ib++) {
    xwork->data[ib] = 0.0;
  }

  nNaNs = 0;
  ib = 0;
  for (k = 0; k + 1 <= n; k++) {
    if (rtIsNaN(b_x->data[k])) {
      idx->data[(n - nNaNs) - 1] = k + 1;
      xwork->data[(n - nNaNs) - 1] = b_x->data[k];
      nNaNs++;
    } else {
      ib++;
      idx4[ib - 1] = k + 1;
      x4[ib - 1] = b_x->data[k];
      if (ib == 4) {
        ib = k - nNaNs;
        if (x4[0] >= x4[1]) {
          m = 1;
          wOffset = 2;
        } else {
          m = 2;
          wOffset = 1;
        }

        if (x4[2] >= x4[3]) {
          p = 3;
          i4 = 4;
        } else {
          p = 4;
          i4 = 3;
        }

        if (x4[m - 1] >= x4[p - 1]) {
          if (x4[wOffset - 1] >= x4[p - 1]) {
            perm[0] = (signed char)m;
            perm[1] = (signed char)wOffset;
            perm[2] = (signed char)p;
            perm[3] = (signed char)i4;
          } else if (x4[wOffset - 1] >= x4[i4 - 1]) {
            perm[0] = (signed char)m;
            perm[1] = (signed char)p;
            perm[2] = (signed char)wOffset;
            perm[3] = (signed char)i4;
          } else {
            perm[0] = (signed char)m;
            perm[1] = (signed char)p;
            perm[2] = (signed char)i4;
            perm[3] = (signed char)wOffset;
          }
        } else if (x4[m - 1] >= x4[i4 - 1]) {
          if (x4[wOffset - 1] >= x4[i4 - 1]) {
            perm[0] = (signed char)p;
            perm[1] = (signed char)m;
            perm[2] = (signed char)wOffset;
            perm[3] = (signed char)i4;
          } else {
            perm[0] = (signed char)p;
            perm[1] = (signed char)m;
            perm[2] = (signed char)i4;
            perm[3] = (signed char)wOffset;
          }
        } else {
          perm[0] = (signed char)p;
          perm[1] = (signed char)i4;
          perm[2] = (signed char)m;
          perm[3] = (signed char)wOffset;
        }

        idx->data[ib - 3] = idx4[perm[0] - 1];
        idx->data[ib - 2] = idx4[perm[1] - 1];
        idx->data[ib - 1] = idx4[perm[2] - 1];
        idx->data[ib] = idx4[perm[3] - 1];
        b_x->data[ib - 3] = x4[perm[0] - 1];
        b_x->data[ib - 2] = x4[perm[1] - 1];
        b_x->data[ib - 1] = x4[perm[2] - 1];
        b_x->data[ib] = x4[perm[3] - 1];
        ib = 0;
      }
    }
  }

  wOffset = (x->size[0] - nNaNs) - 1;
  if (ib > 0) {
    for (m = 0; m < 4; m++) {
      perm[m] = 0;
    }

    if (ib == 1) {
      perm[0] = 1;
    } else if (ib == 2) {
      if (x4[0] >= x4[1]) {
        perm[0] = 1;
        perm[1] = 2;
      } else {
        perm[0] = 2;
        perm[1] = 1;
      }
    } else if (x4[0] >= x4[1]) {
      if (x4[1] >= x4[2]) {
        perm[0] = 1;
        perm[1] = 2;
        perm[2] = 3;
      } else if (x4[0] >= x4[2]) {
        perm[0] = 1;
        perm[1] = 3;
        perm[2] = 2;
      } else {
        perm[0] = 3;
        perm[1] = 1;
        perm[2] = 2;
      }
    } else if (x4[0] >= x4[2]) {
      perm[0] = 2;
      perm[1] = 1;
      perm[2] = 3;
    } else if (x4[1] >= x4[2]) {
      perm[0] = 2;
      perm[1] = 3;
      perm[2] = 1;
    } else {
      perm[0] = 3;
      perm[1] = 2;
      perm[2] = 1;
    }

    for (k = 1; k <= ib; k++) {
      idx->data[(wOffset - ib) + k] = idx4[perm[k - 1] - 1];
      b_x->data[(wOffset - ib) + k] = x4[perm[k - 1] - 1];
    }
  }

  m = nNaNs >> 1;
  for (k = 1; k <= m; k++) {
    ib = idx->data[wOffset + k];
    idx->data[wOffset + k] = idx->data[n - k];
    idx->data[n - k] = ib;
    b_x->data[wOffset + k] = xwork->data[n - k];
    b_x->data[n - k] = xwork->data[wOffset + k];
  }

  if ((nNaNs & 1) != 0) {
    b_x->data[(wOffset + m) + 1] = xwork->data[(wOffset + m) + 1];
  }

  nNonNaN = x->size[0] - nNaNs;
  m = 2;
  if (nNonNaN > 1) {
    if (x->size[0] >= 256) {
      nBlocks = nNonNaN >> 8;
      if (nBlocks > 0) {
        for (i4 = 1; i4 <= nBlocks; i4++) {
          n = (i4 - 1) << 8;
          for (b = 0; b < 6; b++) {
            bLen = 1 << (b + 2);
            bLen2 = bLen << 1;
            nPairs = 256 >> (b + 3);
            for (k = 1; k <= nPairs; k++) {
              m = n + (k - 1) * bLen2;
              for (ib = 0; ib + 1 <= bLen2; ib++) {
                b_iwork[ib] = idx->data[m + ib];
                b_xwork[ib] = b_x->data[m + ib];
              }

              p = 0;
              wOffset = bLen;
              ib = m - 1;
              do {
                exitg1 = 0;
                ib++;
                if (b_xwork[p] >= b_xwork[wOffset]) {
                  idx->data[ib] = b_iwork[p];
                  b_x->data[ib] = b_xwork[p];
                  if (p + 1 < bLen) {
                    p++;
                  } else {
                    exitg1 = 1;
                  }
                } else {
                  idx->data[ib] = b_iwork[wOffset];
                  b_x->data[ib] = b_xwork[wOffset];
                  if (wOffset + 1 < bLen2) {
                    wOffset++;
                  } else {
                    ib = (ib - p) + 1;
                    while (p + 1 <= bLen) {
                      idx->data[ib + p] = b_iwork[p];
                      b_x->data[ib + p] = b_xwork[p];
                      p++;
                    }

                    exitg1 = 1;
                  }
                }
              } while (exitg1 == 0);
            }
          }
        }

        m = nBlocks << 8;
        ib = nNonNaN - m;
        if (ib > 0) {
          merge_block(idx, b_x, m, ib, 2, iwork, xwork);
        }

        m = 8;
      }
    }

    merge_block(idx, b_x, 0, nNonNaN, m, iwork, xwork);
  }

  if ((nNaNs > 0) && (nNonNaN > 0)) {
    for (k = 0; k + 1 <= nNaNs; k++) {
      xwork->data[k] = b_x->data[nNonNaN + k];
      iwork->data[k] = idx->data[nNonNaN + k];
    }

    for (k = nNonNaN - 1; k + 1 > 0; k--) {
      b_x->data[nNaNs + k] = b_x->data[k];
      idx->data[nNaNs + k] = idx->data[k];
    }

    for (k = 0; k + 1 <= nNaNs; k++) {
      b_x->data[k] = xwork->data[k];
      idx->data[k] = iwork->data[k];
    }
  }

  emxFree_real_T(&xwork);
  emxFree_int32_T(&iwork);
  ib = x->size[0];
  x->size[0] = b_x->size[0];
  emxEnsureCapacity((emxArray__common *)x, ib, (int)sizeof(double));
  m = b_x->size[0];
  for (ib = 0; ib < m; ib++) {
    x->data[ib] = b_x->data[ib];
  }

  emxFree_real_T(&b_x);
}

//
// Arguments    : const double x[50]
// Return Type  : double
//
static double trapz(const double x[50])
{
  double z;
  int iy;
  double ylast;
  int k;
  z = 0.0;
  iy = 0;
  ylast = x[0];
  for (k = 0; k < 49; k++) {
    iy++;
    z += (ylast + x[iy]) / 2.0;
    ylast = x[iy];
  }

  return z;
}

//
// Arguments    : emxArray_real_T *A
//                emxArray_real_T *tau
//                emxArray_int32_T *jpvt
// Return Type  : void
//
static void xgeqp3(emxArray_real_T *A, emxArray_real_T *tau, emxArray_int32_T
                   *jpvt)
{
  int n;
  int mn;
  int k;
  int b_n;
  int yk;
  emxArray_real_T *work;
  emxArray_real_T *vn1;
  emxArray_real_T *vn2;
  int nmi;
  int i;
  emxArray_real_T *x;
  int i_i;
  int ix;
  double smax;
  double atmp;
  double s;
  int i_ip1;
  int lastc;
  boolean_T exitg2;
  int exitg1;
  n = A->size[1];
  if (1 <= A->size[1]) {
    mn = 1;
  } else {
    mn = A->size[1];
  }

  k = tau->size[0];
  tau->size[0] = mn;
  emxEnsureCapacity((emxArray__common *)tau, k, (int)sizeof(double));
  if (A->size[1] < 1) {
    b_n = 0;
  } else {
    b_n = A->size[1];
  }

  k = jpvt->size[0] * jpvt->size[1];
  jpvt->size[0] = 1;
  jpvt->size[1] = b_n;
  emxEnsureCapacity((emxArray__common *)jpvt, k, (int)sizeof(int));
  if (b_n > 0) {
    jpvt->data[0] = 1;
    yk = 1;
    for (k = 2; k <= b_n; k++) {
      yk++;
      jpvt->data[k - 1] = yk;
    }
  }

  if (A->size[1] == 0) {
  } else {
    emxInit_real_T1(&work, 1);
    yk = A->size[1];
    k = work->size[0];
    work->size[0] = yk;
    emxEnsureCapacity((emxArray__common *)work, k, (int)sizeof(double));
    for (k = 0; k < yk; k++) {
      work->data[k] = 0.0;
    }

    emxInit_real_T1(&vn1, 1);
    emxInit_real_T1(&vn2, 1);
    yk = A->size[1];
    k = vn1->size[0];
    vn1->size[0] = yk;
    emxEnsureCapacity((emxArray__common *)vn1, k, (int)sizeof(double));
    k = vn2->size[0];
    vn2->size[0] = vn1->size[0];
    emxEnsureCapacity((emxArray__common *)vn2, k, (int)sizeof(double));
    k = 0;
    for (nmi = 0; nmi + 1 <= n; nmi++) {
      vn1->data[nmi] = std::abs(A->data[k]);
      vn2->data[nmi] = vn1->data[nmi];
      k++;
    }

    i = 0;
    emxInit_real_T(&x, 2);
    while (i + 1 <= mn) {
      i_i = i + i;
      nmi = n - i;
      if (nmi < 1) {
        yk = 0;
      } else {
        yk = 1;
        if (nmi > 1) {
          ix = i;
          smax = vn1->data[i];
          for (k = 2; k <= nmi; k++) {
            ix++;
            s = vn1->data[ix];
            if (s > smax) {
              yk = k;
              smax = s;
            }
          }
        }
      }

      b_n = (i + yk) - 1;
      if (b_n + 1 != i + 1) {
        k = x->size[0] * x->size[1];
        x->size[0] = 1;
        x->size[1] = A->size[1];
        emxEnsureCapacity((emxArray__common *)x, k, (int)sizeof(double));
        yk = A->size[0] * A->size[1];
        for (k = 0; k < yk; k++) {
          x->data[k] = A->data[k];
        }

        x->data[b_n] = A->data[i];
        x->data[i] = A->data[b_n];
        k = A->size[0] * A->size[1];
        A->size[0] = 1;
        A->size[1] = x->size[1];
        emxEnsureCapacity((emxArray__common *)A, k, (int)sizeof(double));
        yk = x->size[1];
        for (k = 0; k < yk; k++) {
          A->data[A->size[0] * k] = x->data[x->size[0] * k];
        }

        yk = jpvt->data[b_n];
        jpvt->data[b_n] = jpvt->data[i];
        jpvt->data[i] = yk;
        vn1->data[b_n] = vn1->data[i];
        vn2->data[b_n] = vn2->data[i];
      }

      if (i + 1 < 1) {
        atmp = A->data[i_i];
        s = 0.0;
        if (1 - i <= 0) {
        } else {
          smax = xnrm2(-i, A, i_i + 2);
          if (smax != 0.0) {
            smax = rt_hypotd_snf(A->data[i_i], smax);
            if (A->data[i_i] >= 0.0) {
              smax = -smax;
            }

            if (std::abs(smax) < 1.0020841800044864E-292) {
              yk = 0;
              do {
                yk++;
                xscal(-i, 9.9792015476736E+291, A, i_i + 2);
                smax *= 9.9792015476736E+291;
                atmp *= 9.9792015476736E+291;
              } while (!(std::abs(smax) >= 1.0020841800044864E-292));

              smax = xnrm2(-i, A, i_i + 2);
              smax = rt_hypotd_snf(atmp, smax);
              if (atmp >= 0.0) {
                smax = -smax;
              }

              s = (smax - atmp) / smax;
              xscal(-i, 1.0 / (atmp - smax), A, i_i + 2);
              for (k = 1; k <= yk; k++) {
                smax *= 1.0020841800044864E-292;
              }

              atmp = smax;
            } else {
              s = (smax - A->data[i_i]) / smax;
              xscal(-i, 1.0 / (A->data[i_i] - smax), A, i_i + 2);
              atmp = smax;
            }
          }
        }

        tau->data[i] = s;
        A->data[i_i] = atmp;
      } else {
        tau->data[i] = 0.0;
      }

      if (i + 1 < n) {
        atmp = A->data[i_i];
        A->data[i_i] = 1.0;
        i_ip1 = (i + i) + 2;
        if (tau->data[i] != 0.0) {
          b_n = 1 - i;
          yk = i_i - i;
          while ((b_n > 0) && (A->data[yk] == 0.0)) {
            b_n = 0;
            yk--;
          }

          lastc = nmi - 1;
          exitg2 = false;
          while ((!exitg2) && (lastc > 0)) {
            yk = (i_ip1 + lastc) - 1;
            nmi = yk;
            do {
              exitg1 = 0;
              if (nmi <= (yk + b_n) - 1) {
                if (A->data[nmi - 1] != 0.0) {
                  exitg1 = 1;
                } else {
                  nmi++;
                }
              } else {
                lastc--;
                exitg1 = 2;
              }
            } while (exitg1 == 0);

            if (exitg1 == 1) {
              exitg2 = true;
            }
          }
        } else {
          b_n = 0;
          lastc = 0;
        }

        if (b_n > 0) {
          if (lastc == 0) {
          } else {
            for (yk = 1; yk <= lastc; yk++) {
              work->data[yk - 1] = 0.0;
            }

            yk = 0;
            k = (i_ip1 + lastc) - 1;
            for (b_n = i_ip1; b_n <= k; b_n++) {
              ix = i_i;
              smax = 0.0;
              for (nmi = b_n; nmi <= b_n; nmi++) {
                smax += A->data[nmi - 1] * A->data[ix];
                ix++;
              }

              work->data[yk] += smax;
              yk++;
            }
          }

          if (-tau->data[i] == 0.0) {
          } else {
            yk = 0;
            for (nmi = 1; nmi <= lastc; nmi++) {
              if (work->data[yk] != 0.0) {
                smax = work->data[yk] * -tau->data[i];
                ix = i_i;
                for (b_n = i_ip1; b_n <= i_ip1; b_n++) {
                  A->data[b_n - 1] += A->data[ix] * smax;
                  ix++;
                }
              }

              yk++;
              i_ip1++;
            }
          }
        }

        A->data[i_i] = atmp;
      }

      for (nmi = i + 1; nmi + 1 <= n; nmi++) {
        if (vn1->data[nmi] != 0.0) {
          smax = std::abs(A->data[A->size[0] * nmi]) / vn1->data[nmi];
          smax = 1.0 - smax * smax;
          if (smax < 0.0) {
            smax = 0.0;
          }

          s = vn1->data[nmi] / vn2->data[nmi];
          s = smax * (s * s);
          if (s <= 1.4901161193847656E-8) {
            vn1->data[nmi] = 0.0;
            vn2->data[nmi] = 0.0;
          } else {
            vn1->data[nmi] *= std::sqrt(smax);
          }
        }
      }

      i++;
    }

    emxFree_real_T(&x);
    emxFree_real_T(&vn2);
    emxFree_real_T(&vn1);
    emxFree_real_T(&work);
  }
}

//
// Arguments    : int n
//                const emxArray_real_T *x
//                int ix0
// Return Type  : double
//
static double xnrm2(int n, const emxArray_real_T *x, int ix0)
{
  double y;
  double scale;
  int kend;
  int k;
  double absxk;
  double t;
  y = 0.0;
  if (n < 1) {
  } else if (n == 1) {
    y = std::abs(x->data[ix0 - 1]);
  } else {
    scale = 2.2250738585072014E-308;
    kend = (ix0 + n) - 1;
    for (k = ix0; k <= kend; k++) {
      absxk = std::abs(x->data[k - 1]);
      if (absxk > scale) {
        t = scale / absxk;
        y = 1.0 + y * t * t;
        scale = absxk;
      } else {
        t = absxk / scale;
        y += t * t;
      }
    }

    y = scale * std::sqrt(y);
  }

  return y;
}

//
// Arguments    : int n
//                double a
//                emxArray_real_T *x
//                int ix0
// Return Type  : void
//
static void xscal(int n, double a, emxArray_real_T *x, int ix0)
{
  int i9;
  int k;
  i9 = (ix0 + n) - 1;
  for (k = ix0; k <= i9; k++) {
    x->data[k - 1] *= a;
  }
}

//
// Arguments    : const double accX[50]
//                const double accY[50]
//                const double accZ[50]
//                const double gyrX[50]
//                const double gyrY[50]
//                const double gyrZ[50]
// Return Type  : double
//
double activityClassifier(const double accX[50], const double accY[50], const
  double accZ[50], const double gyrX[50], const double gyrY[50], const double
  gyrZ[50])
{
  double Y;
  double A[150];
  double G[150];
  int k;
  double AR[50];
  double GR[50];
  boolean_T b0[50];
  boolean_T b1[50];
  double y;
  emxArray_real_T *FAX;
  emxArray_real_T *FAY;
  emxArray_real_T *FAZ;
  emxArray_real_T *FAR;
  emxArray_real_T *FGX;
  emxArray_real_T *FGY;
  emxArray_real_T *FGZ;
  emxArray_real_T *FGR;
  emxArray_real_T *b_FAX;
  int loop_ub;
  Y = 0.0;

  // %%Divide to make into G instead of m/s^2
  //  Only keep last 50 samples!
  // Convert to degrees/s
  for (k = 0; k < 50; k++) {
    A[k] = accX[k] / 9.806;
    A[50 + k] = accY[k] / 9.806;
    A[100 + k] = accZ[k] / 9.806;
    G[k] = gyrX[k] * 57.2958;
    G[50 + k] = gyrY[k] * 57.2958;
    G[100 + k] = gyrZ[k] * 57.2958;
  }

  resultant(A, AR);
  resultant(G, GR);

  // %%%%%%%%%%%%%%%%%%%%%%%%
  //         .::::::.       %
  //       .::::::::::.     %
  //     .: :::::::::: :.   %
  //  .:.:: THRESHOLDS ::.:.%
  // %%%%%%%%%%%%%%%%%%%%%%%%
  // Threshold1
  // Threshold2
  for (k = 0; k < 50; k++) {
    b0[k] = (AR[k] - 1.0 > 0.4);
    b1[k] = (GR[k] > 63.0);
    AR[k]--;
  }

  y = b0[0];
  for (k = 0; k < 49; k++) {
    y += (double)b0[k + 1];
  }

  emxInit_real_T(&FAX, 2);
  emxInit_real_T(&FAY, 2);
  emxInit_real_T(&FAZ, 2);
  emxInit_real_T(&FAR, 2);
  emxInit_real_T(&FGX, 2);
  emxInit_real_T(&FGY, 2);
  emxInit_real_T(&FGZ, 2);
  emxInit_real_T(&FGR, 2);
  emxInit_real_T(&b_FAX, 2);
  if (y != 0.0) {
    y = b1[0];
    for (k = 0; k < 49; k++) {
      y += (double)b1[k + 1];
    }

    if (y != 0.0) {
      //  Training Data :: TODO:
      // Extract Features and Classify
      activityFeatureExtractionAcc(*(double (*)[50])&A[0], FAX);
      activityFeatureExtractionAcc(*(double (*)[50])&A[50], FAY);
      activityFeatureExtractionAcc(*(double (*)[50])&A[100], FAZ);
      activityFeatureExtractionAcc(AR, FAR);
      activityFeatureExtractionGyro(*(double (*)[50])&G[0], FGX);
      activityFeatureExtractionGyro(*(double (*)[50])&G[50], FGY);
      activityFeatureExtractionGyro(*(double (*)[50])&G[100], FGZ);
      activityFeatureExtractionGyro(GR, FGR);

      // Apply some threshold? else ? 0
      k = b_FAX->size[0] * b_FAX->size[1];
      b_FAX->size[0] = 1;
      b_FAX->size[1] = ((((((FAX->size[1] + FAY->size[1]) + FAZ->size[1]) +
                           FAR->size[1]) + FGX->size[1]) + FGY->size[1]) +
                        FGZ->size[1]) + FGR->size[1];
      emxEnsureCapacity((emxArray__common *)b_FAX, k, (int)sizeof(double));
      loop_ub = FAX->size[1];
      for (k = 0; k < loop_ub; k++) {
        b_FAX->data[b_FAX->size[0] * k] = FAX->data[FAX->size[0] * k];
      }

      loop_ub = FAY->size[1];
      for (k = 0; k < loop_ub; k++) {
        b_FAX->data[b_FAX->size[0] * (k + FAX->size[1])] = FAY->data[FAY->size[0]
          * k];
      }

      loop_ub = FAZ->size[1];
      for (k = 0; k < loop_ub; k++) {
        b_FAX->data[b_FAX->size[0] * ((k + FAX->size[1]) + FAY->size[1])] =
          FAZ->data[FAZ->size[0] * k];
      }

      loop_ub = FAR->size[1];
      for (k = 0; k < loop_ub; k++) {
        b_FAX->data[b_FAX->size[0] * (((k + FAX->size[1]) + FAY->size[1]) +
          FAZ->size[1])] = FAR->data[FAR->size[0] * k];
      }

      loop_ub = FGX->size[1];
      for (k = 0; k < loop_ub; k++) {
        b_FAX->data[b_FAX->size[0] * ((((k + FAX->size[1]) + FAY->size[1]) +
          FAZ->size[1]) + FAR->size[1])] = FGX->data[FGX->size[0] * k];
      }

      loop_ub = FGY->size[1];
      for (k = 0; k < loop_ub; k++) {
        b_FAX->data[b_FAX->size[0] * (((((k + FAX->size[1]) + FAY->size[1]) +
          FAZ->size[1]) + FAR->size[1]) + FGX->size[1])] = FGY->data[FGY->size[0]
          * k];
      }

      loop_ub = FGZ->size[1];
      for (k = 0; k < loop_ub; k++) {
        b_FAX->data[b_FAX->size[0] * ((((((k + FAX->size[1]) + FAY->size[1]) +
          FAZ->size[1]) + FAR->size[1]) + FGX->size[1]) + FGY->size[1])] =
          FGZ->data[FGZ->size[0] * k];
      }

      loop_ub = FGR->size[1];
      for (k = 0; k < loop_ub; k++) {
        b_FAX->data[b_FAX->size[0] * (((((((k + FAX->size[1]) + FAY->size[1]) +
          FAZ->size[1]) + FAR->size[1]) + FGX->size[1]) + FGY->size[1]) +
          FGZ->size[1])] = FGR->data[FGR->size[0] * k];
      }

      Y = knnclassify(b_FAX);
    }
  }

  emxFree_real_T(&b_FAX);
  emxFree_real_T(&FGR);
  emxFree_real_T(&FGZ);
  emxFree_real_T(&FGY);
  emxFree_real_T(&FGX);
  emxFree_real_T(&FAR);
  emxFree_real_T(&FAZ);
  emxFree_real_T(&FAY);
  emxFree_real_T(&FAX);

  // - function activityClassifier -%
  return Y;
}

//
// Arguments    : void
// Return Type  : void
//
void activityClassifier_initialize()
{
  rt_InitInfAndNaN(8U);
}

//
// Arguments    : void
// Return Type  : void
//
void activityClassifier_terminate()
{
}

//
// File trailer for activityClassifier.cpp
//
// [EOF]
//
