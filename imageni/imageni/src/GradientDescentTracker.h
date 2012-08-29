/*
 * Tracker.h
 *
 *  Created on: Jul 22, 2012
 *      Author: joohwile
 */

#ifndef TRACKER_H_
#define TRACKER_H_

#include "config.h"
#include "MathCode.h"
#include "Image.h"
#include "LinearTransform.h"
#include "Gradient.h"
#include "ProgressReport.h"
#include "Clock.h"

#define PARAM_DIMENSION 4
template<typename Real, typename Pixel>
class CGradientDescentTracker {
public:
	typedef MatrixR<double> MatrixType;
	typedef CImage<Pixel> ImageType;
	typedef CImage<Real> RealImageType;
private:
	ProgressReport* _report;
	int _iter;
	ImageType* _source;
	RealImageType* _target;
	RealImageType* _residual;
	Real _threshold;
	RealImageType _I1, _I2;
	CImageResample _resample;
	Real _cost;
	int _w;
	int _h;
	int _n;
public:
	CGradientDescentTracker(ImageType* t) {
		_source = t;
		_w = _source->getWidth();
		_h = _source->getHeight();
		_n = _w * _h;
		_target = NULL;
		_residual = NULL;
		_iter = 0;
		_report = NULL;
		_cost = 0;
	}
	virtual ~CGradientDescentTracker() {
	}
	void precompute(MatrixType &M0tM0invM0t) {
	}
	void setReporter(ProgressReport* reporter) {
		_report = reporter;
	}
	void setPrior(MatrixType* A) {
	}
	void setCostThreshold(Real threshold) {
		_threshold = threshold;
	}

private:
	Real computeResidual(ImageType* source, RealImageType* target,
			RealImageType* residual) {
		int n = source->getWidth() * source->getHeight();
		Pixel* src1 = source->getData();
		Real* src2 = target->getData();
		Real* dst = residual->getData();
		Real cost = 0;
		for (int i = 0; i < n; i++) {
		    if (abs(src2[i]) > 255 || abs(src1[i]) > 255) {
		        LOGD("Wrong data : %f, %d", src2[i], src1[i]);
		        return 1e9;
		    }
			dst[i] = (src2[i] - src1[i]);
			cost += (dst[i] * dst[i]);
		}
		cost /= n;
		if (cost > 255*255) {
		    LOGD("Wrong cost! n = %d, (s1, s2, d)[100] = (%d, %f, %f)", n, src1[100], src2[100], dst[100]);
		}
		//LOGD("cost = %f, (src1, src2, dst)[100] = %d, %f, %f", cost, src1[100], src2[100], dst[100]);
		return cost;
	}

	inline void resampleImage(ImageType& frame, TransformParamType& x,
			RealImageType& capture) {
		Mat3 m;
		if (PARAM_DIMENSION == 3) {
			ImageTransform::createRigidTransformForImage(x, _w, _h, m);
		} else if (PARAM_DIMENSION == 4) {
			ImageTransform::createSimilarityTransformForImage(x, _w, _h, m);
		}
		_resample.resampleIJ2XYinReal<Pixel, Real>(frame, *_source, m, capture);
	}

	/**
	 * 1) compute residual and cost
	 * 2) compute jacobian of S(\mu)
	 */
	bool optimize(ImageType& frame, TransformParamType& x0,
			TransformParamType& xs) {
		Vec4d delta;
		float steps[4] = { 1, 1, 1, 1 };
		for (int i = 0; i < PARAM_DIMENSION; i++) {
			TransformParamType x1(x0), x2(x0);
			x1._x[i] += steps[i];
			x2._x[i] -= steps[i];
			resampleImage(frame, x1, _I1);
			resampleImage(frame, x2, _I2);
			Real* I1 = _I1.getData();
			Real* I2 = _I2.getData();
			Real dx = 0;
			for (int j = 0; j < _n; j++) {
				dx += (*_residual)[j] * (I1[j] - I2[j]);
			}
			delta[i] = dx;
		}

		if (delta.squaredSum() == 0) {
			xs.copyFrom(x0);
			return false;
		}

		Vec4d grad = delta.normal();
		grad = grad.round();
		xs.copyFrom(x0);
		for (int i = 0; i < PARAM_DIMENSION; i++) {
			xs._x[i] = x0._x[i] - grad[i];
		}

		double delta2 = delta.squaredSum();
		return delta2 > 0.1;
	}

public:
	Real getCost() {
		return _cost;
	}

	int track(ImageType& frame, TransformParamType x0, TransformParamType& xs,
			int maxIter) {
		RealImageType trackingTarget;
		trackingTarget.createEmptyOf(_w, _h);
		_target = &trackingTarget;

		RealImageType trackingResidual;
		trackingResidual.createEmptyOf(_w, _h);
		_residual = &trackingResidual;

		_I1.createEmptyOf(_w, _h);
		_I2.createEmptyOf(_w, _h);

		bool cont = true;
		Real prevCost = 1e9;
		for (_iter = 0; _iter < maxIter && cont; _iter++) {
			resampleImage(frame, x0, *_target);
			_cost = computeResidual(_source, _target, _residual);
			if (_report != NULL) {
				_report->report(_iter, _cost);
			}
			if (prevCost < _cost) {
				break;
			}
			cont = optimize(frame, x0, xs);
			x0.copyFrom(xs);
			prevCost = _cost;
		}
		xs.copyFrom(x0);
		return 0;
	}
};
#endif /* TRACKER_H_ */

