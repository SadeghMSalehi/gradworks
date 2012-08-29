/*
 * NativeCode.cpp
 *
 *  Created on: Aug 10, 2012
 *      Author: joohwile
 */

#include "NativeCode.hpp"
#include "MathCode.h"
#include "Image.h"
#include "HaarFaceDetector.h"
#include "LinearTransform.h"
#include "jnihelper.h"

typedef CImage<jint> ImageType;

static jobject newTransformParam(JNIEnv* env, TransformParamType params) {
	jclass txfParamCls =
			safeFindClassMacro(txfParamCls, "Lcom/intel/vpg/MotionParameter;")
	jmethodID txfParamClsCons =
			safeFindMethodMacro(txfParamClsCons, txfParamCls, "<init>", "()V")
	jobject xs = env->NewObject(txfParamCls, txfParamClsCons);
	jfieldID txFld = env->GetFieldID(txfParamCls, "tx", "F");
	jfieldID tyFld = env->GetFieldID(txfParamCls, "ty", "F");
	jfieldID thetaFld = env->GetFieldID(txfParamCls, "thetaInDegree", "F");
	jfieldID scaleFld = env->GetFieldID(txfParamCls, "scalePercentage", "F");
	env->SetFloatField(xs, txFld, params.tx);
	env->SetFloatField(xs, tyFld, params.ty);
	env->SetFloatField(xs, thetaFld, params.theta);
	env->SetFloatField(xs, scaleFld, params.scale);
	return xs;
}

static bool useTransformParam(JNIEnv* env, jobject jxs,
		TransformParamType& xs) {
	jclass txfParamCls =
			safeFindClassMacro(txfParamCls, "Lcom/intel/vpg/MotionParameter;")
	jfieldID txFld = env->GetFieldID(txfParamCls, "tx", "F");
	jfieldID tyFld = env->GetFieldID(txfParamCls, "ty", "F");
	jfieldID thetaFld = env->GetFieldID(txfParamCls, "thetaInDegree", "F");
	jfieldID scaleFld = env->GetFieldID(txfParamCls, "scalePercentage", "F");

	xs.tx = env->GetFloatField(jxs, txFld);
	xs.ty = env->GetFloatField(jxs, tyFld);
	xs.theta = env->GetFloatField(jxs, thetaFld);
	xs.scale = env->GetFloatField(jxs, scaleFld);
	return true;
}

static jobject newMatrix3(JNIEnv* env, Mat3 params) {
	jclass cls = safeFindClassMacro(cls, "Lcom/intel/vpg/Matrix;")
	jmethodID clsCons = safeFindMethodMacro(clsCons, cls, "<init>", "(II)V")
	jobject mat = env->NewObject(cls, clsCons, 3, 3);
	jfieldID dataField = env->GetFieldID(cls, "data", "[F");
	jfloatArray jdata = (jfloatArray) env->GetObjectField(mat, dataField);
	arrayCopyToJNIMacro(data, Float, jfloat, params._d)
	return mat;
}

static void useMatrix3(JNIEnv* env, jobject jmat, Mat3& mat) {
	jfieldID dataField = env->GetFieldID(env->GetObjectClass(jmat), "data",
			"[F");
	jfloatArray jdata = (jfloatArray) env->GetObjectField(jmat, dataField);
	arrayCopyFromJNIMacro(data, Float, jfloat, mat._d);
}

static jobject newIntImage(JNIEnv* env, ImageType& image) {
	jclass cls = safeFindClassMacro(cls, "Lcom/intel/vpg/IntImage;")
	jmethodID clsCons = safeFindMethodMacro(clsCons, cls, "<init>", "(II)V")
	jobject jintImage = env->NewObject(cls, clsCons, image.getWidth(),
			image.getHeight());
	jfieldID dataField = env->GetFieldID(cls, "data", "[I");
	jintArray jdata = (jintArray) env->GetObjectField(jintImage, dataField);
	arrayCopyToJNIMacro(data, Int, jint, image.getData())
	return jintImage;
}

static jintArray intImageArray(JNIEnv* env, jobject intImage, int *w, int *h) {
	jclass intImageClass = env->GetObjectClass(intImage);
	jfieldID dataField = env->GetFieldID(intImageClass, "data", "[I");
	jfieldID widthField = env->GetFieldID(intImageClass, "_width", "I");
	jfieldID heightField = env->GetFieldID(intImageClass, "_height", "I");

	jintArray jdata = (jintArray) env->GetObjectField(intImage, dataField);
	jint jwidth = env->GetIntField(intImage, widthField);
	jint jheight = env->GetIntField(intImage, heightField);

	*w = jwidth;
	*h = jheight;
	return jdata;
}

/**
 * must release array data after the use of it
 */
static jintArray useImage(JNIEnv* env, jobject jimg, ImageType& img,
		jint** pdata) {
	int w = 0;
	int h = 0;

	jintArray jdata = intImageArray(env, jimg, &w, &h);
	*pdata = lockArray(data, Int);

	img.wrappingOf(*pdata, w, h);
	return jdata;
}

JNIEXPORT jintArray JNICALL Java_com_intel_vpg_NativeCode_detectFaces(
		JNIEnv *env, jclass that, jobject frameImage, jint minScaleFactor, jboolean earlyReturn) {
	int w = 0;
	int h = 0;

	jintArray jframe = intImageArray(env, frameImage, &w, &h);
	jint* pframe = env->GetIntArrayElements(jframe, NULL);
	CImage<jint> frame(pframe, w, h), faceImg;

	jfloat lastScale = 0;
	CHaarFaceDetector<jint, jfloat> detector;
	detector.setMinScaleFactor((int) minScaleFactor);
	detector.detectAtScales(frame, frame, lastScale, earlyReturn);

	haar_point_array_t faces = detector.getResult();
	jintArray jfaces = env->NewIntArray(faces.size() * 4);
	jint* pfaces = lockArray(faces, Int);
	for (size_t i = 0; i < faces.size(); i++) {
		pfaces[4 * i] = faces[i].x;
		pfaces[4 * i + 1] = faces[i].y;
		pfaces[4 * i + 2] = faces[i].w;
		pfaces[4 * i + 3] = faces[i].h;
	}
	releaseArray(faces, Int);
	return jfaces;
}

/*
 * Class:     com_intel_vpg_NativeCode
 * Method:    resampleRegion
 * Signature: (Lcom/intel/vpg/IntImage;Lcom/intel/vpg/IntImage;Lcom/intel/vpg/MotionParameter;)Lcom/intel/vpg/IntImage;
 */JNIEXPORT jobject JNICALL Java_com_intel_vpg_NativeCode_resampleRegion(
		JNIEnv *env, jclass that, jobject jframe, jobject jobj,
		jobject jtxfParam) {

	ImageType frame, obj, dst;

	jint *pframeArray, *pobjArray;
	jintArray jframeArray = useImage(env, jframe, frame, &pframeArray);
	jintArray jobjArray = useImage(env, jobj, obj, &pobjArray);

	Mat3 txf;

	jclass tpCls = env->GetObjectClass(jtxfParam);
	jfieldID txfId = env->GetFieldID(env->GetObjectClass(jtxfParam),
			"forwardTransform", "Lcom/intel/vpg/Matrix;");
	jobject jtxfMat = env->GetObjectField(jtxfParam, txfId);
	useMatrix3(env, jtxfMat, txf);

	CImageResample resample;
	resample.resampleIJ2XYinReal<jint,jint>(frame, obj, txf, dst);

	jobject joutImage = newIntImage(env, dst);
	releaseArray(frameArray, Int);
	releaseArray(objArray, Int);

	return joutImage;
}

/*
 * Class:     com_intel_vpg_NativeCode
 * Method:    resampleRegion
 * Signature: (Lcom/intel/vpg/IntImage;Lcom/intel/vpg/IntImage;Lcom/intel/vpg/MotionParameter;)Lcom/intel/vpg/IntImage;
 */
 JNIEXPORT void JNICALL Java_com_intel_vpg_NativeCode_resize
   (JNIEnv *env, jclass that, jintArray jsrc, jint sw, jint sh, jintArray jdst, jint dw, jint dh) {

	 jint* psrc = lockArray(src, Int);
	 CImage<jint> src(psrc, sw, sh);
	 CImage<jint> dst;
	 dst.createResized(src, dw, dh);
	 arrayCopyToJNIMacro(dst, Int, jint, dst.getData());
	 releaseArray(src, Int);
}
