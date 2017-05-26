//
//

// Include Files
#include "rt_nonfinite.h"
#include "resultant.h"
/*Additional Includes*/
#include <jni.h>
#include <android/log.h>

#define  LOG_TAG "resultant-cpp"
#define  LOGI(...)  __android_log_print(ANDROID_LOG_INFO, LOG_TAG, __VA_ARGS__)
#define  LOGD(...)  __android_log_print(ANDROID_LOG_DEBUG, LOG_TAG, __VA_ARGS__)
#define  LOGE(...)  __android_log_print(ANDROID_LOG_ERROR, LOG_TAG, __VA_ARGS__)

extern "C" {
JNIEXPORT jdouble JNICALL
//Call this function with (data, data, data, data, datalen, Fs);
//Don't need array size; can check size array in C.
Java_musamahmood_fitnessapp_MainActivity_jniMean(JNIEnv *env, jobject jobject1, jdoubleArray data, jint len) {
    jdouble  *dataArray;
    dataArray = env->GetDoubleArrayElements(data, NULL);
    if(dataArray==NULL) LOGE("Error!");
    double result = minDataMean(dataArray, len);
    return result;
}
}
//JNI function call:
extern "C" {
JNIEXPORT jdouble JNICALL
//Call this function with (data, data, data, data, datalen, Fs);
//Don't need array size; can check size array in C.
Java_musamahmood_fitnessapp_MainActivity_jniResultantAcc(JNIEnv *env, jobject jobject1, jdoubleArray data) {
    jdouble  *dataArray;
    dataArray = env->GetDoubleArrayElements(data, NULL);
    if(dataArray==NULL) LOGE("Error!");
    double result = resultantAcc(dataArray);
    return result;
}
}

extern "C" {
JNIEXPORT jdouble JNICALL
//Call this function with (data, data, data, data, datalen, Fs);
//Don't need array size; can check size array in C.
Java_musamahmood_fitnessapp_MainActivity_jniResultantGyro(JNIEnv *env, jobject jobject1, jdoubleArray data) {
    jdouble  *dataArray;
    dataArray = env->GetDoubleArrayElements(data, NULL);
    double result = resultantGyro(dataArray);
    return result;
}
}


// Function Definitions

//Mean 1 minute of data

double minDataMean(const double X[1000], const int len) {
    double X1 = 0;
    for (int i = 0; i < len; ++i) {
        X1+=X[i]; //sum
    }
    X1 /= len; //take average from length
    return X1;
}

//
// Get Resultant. Size of input vector must be exactly 3.
//
// Arguments    : const double X[3]
// Return Type  : double
//
double resultantAcc(const double X[3]) {
  double X1[3];
  int k;
  for (k = 0; k < 3; k++) {
    //CONVERT FROM m/s^2 to G-s
    X1[k] = X[k]/9.806;
    X1[k] = X1[k] * X1[k];
  }

  return std::sqrt((X1[0] + X1[1]) + X1[2]);
}

double resultantGyro(const double X[3]) {
    double X1[3];
    int k;
    for (k = 0; k < 3; k++) {
        X1[k] = X[k] * X[k];
    }
    return std::sqrt((X1[0] + X1[1]) + X1[2]);
}

//
// Arguments    : void
// Return Type  : void
//
void resultant_initialize() {
  rt_InitInfAndNaN(8U);
}


//
// File trailer for resultant.cpp
//
// [EOF]
//
