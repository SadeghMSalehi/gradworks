#include "jni.h"

#include "config.h"

#include <vector>
#include "Image.h"
#include "HaarFaceDetector.h"

typedef unsigned char uint8_t;

typedef struct {
    uint8_t alpha;
    uint8_t red;
    uint8_t green;
    uint8_t blue;
} argb;

extern "C" {
JNIEXPORT void JNICALL Java_com_intel_NativeCall_resizeImage(JNIEnv* env, jclass clazz, jintArray grayArray, jint width,
                jint height, jintArray output, jint newWidth, jint newHeight) {
    jint* inputBuffer = (jint*) env->GetIntArrayElements(grayArray, NULL);
    jint* outputBuffer = (jint*) env->GetIntArrayElements(output, NULL);
    CImage<jint> inputImg(inputBuffer, width, height);
    CImage<jint> outputImg(outputBuffer, newWidth, newHeight);
    outputImg.createResized(inputImg, newWidth, newHeight, 1);
    env->ReleaseIntArrayElements(grayArray, inputBuffer, 0);
    env->ReleaseIntArrayElements(output, outputBuffer, 0);
}


JNIEXPORT void JNICALL Java_com_intel_NativeCall_decodeYUV2Gray(JNIEnv *env, jclass clazz, jbyteArray yuvArray,
                jint width, jint height, jintArray grayArray) {
    jboolean isCopy;
    jsize frameSize = width * height; //env->GetArrayLength(yuvArray);
    jbyte* yuvData = env->GetByteArrayElements(yuvArray, &isCopy);
    jint* grayData = env->GetIntArrayElements(grayArray, &isCopy);

    for (int pix = 0; pix < frameSize; pix++) {
        jint pixVal = (0xff & ((jint) yuvData[pix])); // - 16;
        grayData[pix] = pixVal;
    }

    env->ReleaseByteArrayElements(yuvArray, yuvData, 0);
    env->ReleaseIntArrayElements(grayArray, grayData, 0);
}

JNIEXPORT void JNICALL Java_com_intel_NativeCall_convertGrayToRGB(JNIEnv *env, jclass clazz, jintArray grayArray,
                jint width, jint height, jintArray rgbArray) {
    jboolean isCopy;
    jsize frameSize = width * height; //env->GetArrayLength(yuvArray);
    jint* grayData = env->GetIntArrayElements(grayArray, &isCopy);
    jint* rgbData = env->GetIntArrayElements(rgbArray, &isCopy);

    LOGD("convertGrayToRGB()");
    if (rgbData != NULL) {
        for (int i = 0; i < frameSize; i++) {
            jint argb = 0xff000000 | grayData[i] << 16 | grayData[i] << 8 | grayData[i];
            rgbData[i] = argb;
        }
    }

    env->ReleaseIntArrayElements(grayArray, grayData, 0);
    env->ReleaseIntArrayElements(rgbArray, rgbData, 0);
}

JNIEXPORT void JNICALL Java_com_intel_NativeCall_decodeGray2ARGB(JNIEnv *env, jclass clazz, jbyteArray grayArray,
                jint width, jint height, jintArray argbArray) {
    jboolean isCopy;
    jsize frameSize = width * height;
    jbyte* grayData = env->GetByteArrayElements(grayArray, &isCopy);
    jint* argbData = env->GetIntArrayElements(argbArray, &isCopy);

    for (int i = 0; i < frameSize; i++) {
        jint argb = 0xff000000 | grayData[i] << 16 | grayData[i] << 8 | grayData[i];
        argbData[i] = argb;
    }

    env->ReleaseByteArrayElements(grayArray, grayData, 0);
    env->ReleaseIntArrayElements(argbArray, argbData, 0);

}

JNIEXPORT void JNICALL Java_com_intel_NativeCall_decodeYUV2RGB(JNIEnv *env, jclass clazz, jbyteArray yuvArray,
                jint width, jint height, jintArray argbArray) {
    jboolean isCopy;
    jsize frameSize = width * height; //env->GetArrayLength(yuvArray);
    //LOGI( "(%s:%d) YUV array size: %d (%d x %d)", __FILE__, __LINE__, frameSize, width, height);

    jbyte* yuv420sp = env->GetByteArrayElements(yuvArray, &isCopy);
    //LOGI("yuvData jbyte* acquired");

    jint* argbData = env->GetIntArrayElements(argbArray, &isCopy);
    //LOGI("argbData jint* acquired");

    for (int j = 0, yp = 0; j < height; j++) {
        int uvp = frameSize + (j >> 1) * width, u = 0, v = 0;
        for (int i = 0; i < width; i++, yp++) {
            int y = (0xff & ((int) yuv420sp[yp])) - 16;
            if (y < 0)
                y = 0;
            if ((i & 1) == 0) {
                v = (0xff & yuv420sp[uvp++]) - 128;
                u = (0xff & yuv420sp[uvp++]) - 128;
            }

            int y1192 = 1192 * y;
            int r = (y1192 + 1634 * v);
            int g = (y1192 - 833 * v - 400 * u);
            int b = (y1192 + 2066 * u);

            if (r < 0)
                r = 0;
            else if (r > 262143)
                r = 262143;
            if (g < 0)
                g = 0;
            else if (g > 262143)
                g = 262143;
            if (b < 0)
                b = 0;
            else if (b > 262143)
                b = 262143;

            int rgbValue = 0xff000000 | ((r << 6) & 0xff0000) | ((g >> 2) & 0xff00) | ((b >> 10) & 0xff);
            argbData[yp] = rgbValue;
        }
    }

    env->ReleaseByteArrayElements(yuvArray, yuv420sp, 0);
    env->ReleaseIntArrayElements(argbArray, argbData, 0);
}

}
