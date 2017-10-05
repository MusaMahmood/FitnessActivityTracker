//
// Created by mahmoodms on 4/22/2017.
//
/*Additional Includes for JNI Android Integration*/
#include <jni.h>
#include <android/log.h>
#include "activityClassifier.h"

#define  LOG_TAG "activityTracker-cpp"
#define  LOGI(...)  __android_log_print(ANDROID_LOG_INFO, LOG_TAG, __VA_ARGS__)
#define  LOGD(...)  __android_log_print(ANDROID_LOG_DEBUG, LOG_TAG, __VA_ARGS__)
#define  LOGE(...)  __android_log_print(ANDROID_LOG_ERROR, LOG_TAG, __VA_ARGS__)

// Function Declarations
static void argInit_50x1_real_T(double result[50]);

static double argInit_real_T();

static double main_activityClassifier();

// Function Definitions

//
// Arguments    : double result[50]
// Return Type  : void
//
static void argInit_50x1_real_T(double result[50]) {
    int idx0;

    // Loop over the array to initialize each element.
    for (idx0 = 0; idx0 < 50; idx0++) {
        // Set the value of the array element.
        // Change this value to the value that the application requires.
        result[idx0] = argInit_real_T();
    }
}

//
// Arguments    : void
// Return Type  : double
//
static double argInit_real_T() {
    return 0.0;
}

//
// Arguments    : void
// Return Type  : void
//
static double main_activityClassifier() {
    double dv8[50];
    double dv9[50];
    double dv10[50];
    double dv11[50];
    double dv12[50];
    double dv13[50];
    double Y;

    // Initialize function 'activityClassifier' input arguments.
    // Initialize function input argument 'accX'.
    // Initialize function input argument 'accY'.
    // Initialize function input argument 'accZ'.
    // Initialize function input argument 'gyrX'.
    // Initialize function input argument 'gyrY'.
    // Initialize function input argument 'gyrZ'.
    // Call the entry-point 'activityClassifier'.
    argInit_50x1_real_T(dv8);
    argInit_50x1_real_T(dv9);
    argInit_50x1_real_T(dv10);
    argInit_50x1_real_T(dv11);
    argInit_50x1_real_T(dv12);
    argInit_50x1_real_T(dv13);
    Y = activityClassifier(dv8, dv9, dv10, dv11, dv12, dv13);
    return Y;
}

extern "C" {
JNIEXPORT jdouble JNICALL
Java_musamahmood_fitnessapp_MainActivity_jniActivityTrackerInit(JNIEnv *env, jobject jobject1) {
    activityClassifier_initialize();
    return main_activityClassifier();
}
}



extern "C" {
JNIEXPORT jdouble JNICALL
Java_musamahmood_fitnessapp_MainActivity_jniActivityClassifier(JNIEnv *env, jobject jobject1,
                                                               jdoubleArray accX, jdoubleArray accY,
                                                               jdoubleArray accZ,
                                                               jdoubleArray gyrX, jdoubleArray gyrY,
                                                               jdoubleArray gyrZ) {
    jdouble *accXp = env->GetDoubleArrayElements(accX, NULL);
    jdouble *accYp = env->GetDoubleArrayElements(accY, NULL);
    jdouble *accZp = env->GetDoubleArrayElements(accZ, NULL);
    jdouble *gyrXp = env->GetDoubleArrayElements(gyrX, NULL);
    jdouble *gyrYp = env->GetDoubleArrayElements(gyrY, NULL);
    jdouble *gyrZp = env->GetDoubleArrayElements(gyrZ, NULL);
    double Y = activityClassifier(accXp, accYp, accZp, gyrXp, gyrYp, gyrZp);
    return Y;
}
}

