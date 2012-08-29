/*
 * TrackingContext.cpp
 *
 *  Created on: Aug 10, 2012
 *      Author: joohwile
 */

#include  "jnihelper.h"
#include "TrackingContext.hpp"
#include "Image.h"
#include "GradientDescentTracker.h"
#include "LinearTransform.h"

typedef CGradientDescentTracker<jdouble, jint> TrackerType;

/*
 * Class:     com_intel_vpg_TrackingContext
 * Method:    nativeInit
 * Signature: ([III[D)V
 */JNIEXPORT void JNICALL Java_com_intel_vpg_TrackingContext_nativeInit(
		JNIEnv *env, jobject self, jintArray jtarget, jint w, jint h,
		jdoubleArray jprior) {
	jint* ptarget = lockArray(target, Int);
	CImage<jint> targetImage(ptarget, w, h);

	TrackerType::MatrixType prior(4, w * h);
	TrackerType tracker(&targetImage);
	tracker.precompute(prior);
	arrayCopyToJNIMacro(prior, Double, jdouble, prior._M);

	releaseArray(target, Int);
}

class CostReport: public ProgressReport {
	JNIEnv* _env;
	jobject _self;
	jmethodID _method;

public:
	CostReport(JNIEnv* env, jobject self) {
		_env = env;
		_self = self;
		_method = _env->GetMethodID(_env->GetObjectClass(self), "reportCost",
				"(D)V");
		if (_method == NULL) {
			cout << "Can't find method for report cost" << endl;
			return;
		}
	}
	void report(int iter, double cost) {
		if (_method != NULL) {
			_env->CallVoidMethod(_self, _method, cost);
		}
	}
};


/*
 * Class:     com_intel_vpg_TrackingContext
 * Method:    nativeTrack
 * Signature: ([D[III[III[FI[F)D
 */
JNIEXPORT jdouble JNICALL Java_com_intel_vpg_TrackingContext_nativeTrack(
		JNIEnv *env, jobject self, jdoubleArray jprior, jintArray jframe,
		jint fw, jint fh, jintArray jobject, jint ow, jint oh, jfloatArray jin,
		jint maxIter, jfloatArray jout) {

	jint* pframe = lockArray(frame, Int);
	CImage<jint> frame0(pframe, fw, fh), frame1;
	frame1.createCopyOf(frame0);
	releaseArray(frame, Int);

	jint* pobject = lockArray(object, Int);
	CImage<jint> object0(pobject, ow, oh), object1;
	object1.createCopyOf(object0);
	releaseArray(object, Int);

	jfloat* pin = lockArray(in, Float);
	TransformParamType x0, xs;
	x0.tx = pin[0];
	x0.ty = pin[1];
	x0.theta = pin[2];
	x0.scale = pin[3];
	releaseArray(in, Float);

	CostReport reporter(env, self);
	TrackerType tracker(&object1);
	tracker.track(frame1, x0, xs, maxIter);

	jfloat* pout = lockArray(out, Float);
	pout[0] = xs.tx;
	pout[1] = xs.ty;
	pout[2] = xs.theta;
	pout[3] = xs.scale;
	releaseArray(out, Float);

	return tracker.getCost();
}

/*
 * Class:     com_intel_vpg_TrackingContext
 * Method:    nativeResample
 * Signature: ([III[III)Z
 */JNIEXPORT jboolean JNICALL Java_com_intel_vpg_TrackingContext_nativeResample(
		JNIEnv *env, jobject self, jintArray jframe, jint fw, jint fh,
		jfloatArray jparams, jintArray jcapture, jint cw, jint ch) {
	jint* pframe = lockArray(frame, Int);
	jint* pcapture = lockArray(capture, Int);
	jfloat* pparams = lockArray(params, Float);

	CImage<jint> frame(pframe, fw, fh);
	CImage<jint> capture(pcapture, cw, ch);

	Mat3 fTx;
	TransformParamType x;
	x.tx = pparams[0];
	x.ty = pparams[1];
	x.theta = pparams[2];
	x.scale = pparams[3];
	ImageTransform::createSimilarityTransformForImage(x, cw, ch, fTx);

	CImageResample resample;
	resample.resampleIJ2XYinReal<jint,jint>(frame, capture, fTx, capture);

	releaseArray(frame, Int);
	releaseArray(capture, Int);
	releaseArray(params, Float);

	return true;
}

