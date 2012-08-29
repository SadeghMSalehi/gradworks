/*
 * GaussNewtonTracker.h
 *
 *  Created on: Aug 11, 2012
 *      Author: joohwile
 */

#ifndef GAUSSNEWTONTRACKER_H_
#define GAUSSNEWTONTRACKER_H_

#include "MathCode.h"
#include "Image.h"
#include "Opti.h"
#include "LinearTransform.h"
#include "Gradient.h"
#include "ProgressReport.h"

template<typename Real, typename Pixel>
class CGaussNewtonTracker {
public:
	typedef MatrixR<double> MatrixType;
	typedef CImage<Pixel> ImageType;
	typedef CImage<Real> RealImageType;
private:
	MatrixType* _A;
	Real _prevCost;
	Real _currentCost;
	Real _threshold;
	Real _absoluteThreshold;
	ProgressReport* _report;
	int _iter;
	CImageResample _resample;
	int _w;
	int _h;
	int _n;

	ImageType* _I0;
	RealImageType _R;
	RealImageType _I1;

public:
	CGaussNewtonTracker(ImageType* t) :
			_I0(t) {
		_A = NULL;
		_currentCost = 1e9;
		_threshold = 1000;
		_absoluteThreshold = 500;
		_report = NULL;
		_iter = 0;
		_w = _I0->getWidth();
		_h = _I0->getHeight();
		_n = _w * _h;
		_I1.createEmptyOf(_w, _h);
		_R.createEmptyOf(_w, _h);
	}
	void setReporter(ProgressReport* reporter) {
		_report = reporter;
	}
	void setPrior(MatrixType* A) {
		_A = A;
		_currentCost = 1e9;
		_threshold = 500;
	}
	void setCostThreshold(Real threshold) {
		_threshold = threshold;
	}

	void resampleImage(ImageType& frame, TransformParamType& x) {
		Mat3 txf;
		ImageTransform::createSimilarityTransformForImage(x, _w, _h, txf);
		_resample.resampleIJ2XYinReal<Pixel, Real>(frame, (*_I0), txf, _I1);
	}

	void computeResidual() {
		for (int i = 0; i < _n; i++) {
			_R[i] = _I1[i] - (*_I0)[i];
		}
	}

	bool optimize(ImageType& frame, TransformParamType &x0,
			TransformParamType& xs) {

		Mat4d sigmaT, invSigmaT;
		sigmaT.identity();
		float cs = cosf(-x0.thetaRad());
		float si = sinf(-x0.thetaRad());
		sigmaT[0][0] = 100.f / x0.scale * cs;
		sigmaT[0][1] = 100.f / x0.scale * -si;
		sigmaT[1][0] = 100.f / x0.scale * si;
		sigmaT[1][1] = 100.f / x0.scale * cs;
		sigmaT[3][3] = 100.f / x0.scale;
		sigmaT.inverse(invSigmaT);

		Vec4d M0tE;
		M0tE.zero();
		for (int i = 0; i < _A->nCols; i++) {
			M0tE[0] += _A->_R[0][i] * _R[i];
			M0tE[1] += _A->_R[1][i] * _R[i];
			M0tE[2] += _A->_R[2][i] * _R[i];
			M0tE[3] += _A->_R[3][i] * _R[i];
		}

		Vec4d dMu;
		mult(invSigmaT, M0tE, dMu);

		dMu[2] = dMu[2] * 180 / PI;
		dMu[3] = dMu[3] * 100;

		Vec4d grad = dMu;

		for (int i = 0; i < 4; i++) {
			xs._x[i] = x0._x[i] - grad[i];
		}

		cout << "dMu = " << dMu << "; g = " << grad << "; x0 = " << x0 << "; xs = " << xs << endl;

		double deltaMag = dMu.squaredSum();
		if (deltaMag < 0.1) {
			return false;
		}

		return true;
	}

	int track(ImageType& frame, TransformParamType x0, TransformParamType& xs,
			int maxIter) {
		CClock clock;
		clock.tick();

		bool cont = true;
		_iter = 0;
		TransformParamType px;
		px.copyFrom(x0);
		// if optimization fail at the first step return x0 as xs
		xs.copyFrom(x0);

		for (_iter = 0; _iter < maxIter; _iter++) {
			resampleImage(frame, x0);
			computeResidual();
			cont = optimize(frame, x0, xs);
			if (!cont) {
				xs.copyFrom(px);
				break;
			}
			px.copyFrom(x0);
			x0.copyFrom(xs);
		}
		return clock.tock();
	}

	void precompute(MatrixType &M0tM0invM0t) {
		const int w = _I0->getWidth();
		const int h = _I0->getHeight();

		DoubleImage gx, gy;
		gx.createEmptyOf(w, h);
		gy.createEmptyOf(w, h);

		Gradient::gradientX(*_I0, gx);
		Gradient::gradientY(*_I0, gy);

		int hw = _I0->getWidth() / 2;
		int hh = _I0->getHeight() / 2;

		MatrixType M0t(M0tM0invM0t.nRows, M0tM0invM0t.nCols);
		MatrixType coeff(3, 6);

		for (int j = 0; j < h; j++) {
			for (int i = 0; i < w; i++) {
				const float x = i - hw, y = j - hh;
				const int jwi = j * w + i;

				double* jGx = gx.getScanline(j);
				double* jGy = gy.getScanline(j);

				M0t._R[0][jwi] = jGx[i];
				M0t._R[1][jwi] = jGy[i];
				M0t._R[2][jwi] = -y * jGx[i] + x * jGy[i];
				M0t._R[3][jwi] = x * jGx[i] + y * jGy[i];
				/*
				 #ifdef DEBUG
				 if (isnan(jGx[i])) {
				 cout << "NaN detected at jGx[" << i << "," << j << "] = "
				 << jGx[i] << "/" << gx(i, j) << endl;
				 }
				 if (isnan(jGy[i])) {
				 cout << "NaN detected at jGy[" << i << "," << j << "] = "
				 << jGy[i] << endl;
				 }
				 #endif
				 */
			}
		}

		SquareMatrixR<Real, 4> MtM, invMtM;
		MtM.zero();
		for (int i = 0; i < M0t.nRows; i++) {
			for (int j = 0; j < M0t.nRows; j++) {
				for (int k = 0; k < M0t.nCols; k++) {
					MtM[i][j] += M0t._R[i][k] * M0t._R[j][k];
				}
			}
		}
		MtM.inverse(invMtM);

		for (int i = 0; i < M0tM0invM0t.nRows; i++) {
			for (int j = 0; j < M0tM0invM0t.nCols; j++) {
				M0tM0invM0t[i][j] = 0;
				for (int k = 0; k < M0tM0invM0t.nRows; k++) {
					M0tM0invM0t[i][j] += (invMtM._r[i][k] * M0t[k][j]);
				}
			}
		}
	}
};

#endif /* GAUSSNEWTONTRACKER_H_ */
