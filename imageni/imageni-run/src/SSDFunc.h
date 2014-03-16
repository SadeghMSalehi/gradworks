#include "MathCode.h"
#include "LinearTransform.h"


template<typename I, typename V, typename G, typename H>
class CFunc {
public:
    typedef I InputType;
    typedef V ValueOutputType;
    typedef G GradientOutputType;
    typedef H HessianOutputType;
    virtual ~CFunc() {
    }
    virtual void position(InputType& x) = 0;
    virtual void value(ValueOutputType& value) = 0;
    virtual void gradient(GradientOutputType& grad) = 0;
    virtual void hessian(HessianOutputType& hessian) = 0;
};



template<typename T>
class CSSD: public CFunc<TransformParamType, float, TransformParamType,
		SquareMatrixR<float, 8> > {
private:
	SimilarityTransform _txf;
	//Isometry _txf;
	TransformParamType _p;
	Mat3 _T;
	CImage<T>* _I1;
	CImage<T>* _I0;
	typedef CImage<T> Image;
	int _iRange[10];
	int _jRange[10];
	int _kRange[2];
	int _iStep;
	int _jStep;

public:

	CSSD(Image* fixed, Image* moving) {
		_I1 = fixed;
		_I0 = moving;
		_T.identity();

		_iStep = 1;
		_jStep = 1;
		_kRange[0] = 0;
		_kRange[1] = 1;

		/*
		_iRange[0] = 0;
		_iRange[1] = _I0->getWidth();
		_jRange[0] = 0;
		_jRange[1] = _I0->getHeight();
		*/
		_kRange[0] = 0;
		_kRange[1] = 5;
		_iRange[0] = _I0->getWidth() / 3;
		_iRange[1] = _I0->getWidth() * 2 / 3;
		_iRange[2] = 1;
		_iRange[3] = _I0->getWidth() - 1;
		_iRange[4] = _I0->getWidth() / 3;
		_iRange[5] = _I0->getWidth() * 2 / 3;
		_jRange[0] = 1;
		_jRange[1] = _I0->getHeight() / 3;
		_jRange[2] = _I0->getHeight() / 3;
		_jRange[3] = _I0->getHeight() * 2 / 3;
		_jRange[4] = _I0->getHeight() * 2 / 3;
		_jRange[5] = _I0->getHeight() - 1;
	}
	virtual ~CSSD() {
	}
	Mat3* getTransform() {
		return &_T;
	}
	inline void createTransform(TransformParamType& p, Mat3& T) {
		_txf.createPointTransform(p, T);
	}
	void position(TransformParamType& p) {
		_p = p;
		//LinearTransform2D::createIsometricPointTransform(p, _T);
		createTransform(p, _T);
	}
	void value(float& fv) {
		// assume fixed image is much larger
		Image& I1 = (*_I1);
		Image& I0 = (*_I0);
		double ssd = 0;
		for (int k = _kRange[0]; k <= _kRange[1]; k += 2) {
			for (int j = _jRange[k]; j < _jRange[k + 1]; j += _jStep) {
				for (int i = _iRange[k]; i < _iRange[k + 1]; i += _iStep) {
					Vec2 ij(i, j), xy;
					mult(_T, ij, xy);
					//_DBG_(if (i == 10) { cout << ij << "=>" << xy << endl;})
					if (_I1->isInside(xy[0], xy[1])) {
						double f = I0.bilinear(ij[0], ij[1])
								- I1.bilinear(xy[0], xy[1]);
						double f2 = f * f;
						ssd += f2;
					}
				}
			}
		}
		fv = ssd;
	}
	void gradient(TransformParamType& grad) {
		// assume fixed image is much larger
		Image& I1 = (*_I1);
		Image& I0 = (*_I0);

		Mat3 R1, R0;
		TransformParamType dFdr1 = _p;
		TransformParamType dFdr0 = _p;
		dFdr1.theta += 1;
		dFdr0.theta -= 1;

		Mat3 S1, S0;
		TransformParamType dFds1 = _p;
		TransformParamType dFds0 = _p;
		dFds1.scale += 1;
		dFds0.scale -= 1;

		_txf.createPointTransform(dFdr1, R1);
		_txf.createPointTransform(dFdr0, R0);
		_txf.createPointTransform(dFds1, S1);
		_txf.createPointTransform(dFds0, S0);

		float dtx = 0, dty = 0, dAngle = 0, dScale = 0;
		int validCount = 0;

		for (int k = _kRange[0]; k <= _kRange[1]; k += 2) {
			for (int j = _jRange[k]; j < _jRange[k + 1]; j += _jStep) {
				for (int i = _iRange[k]; i < _iRange[k + 1]; i += _iStep) {
					Vec2 ij(i, j), fp, fr1, fr0, fs1, fs0;
					mult(_T, ij, fp);
					mult(R1, ij, fr1);
					mult(R0, ij, fr0);
					mult(S1, ij, fs1);
					mult(S0, ij, fs0);

					//_DBG_(if (i == 10) { cout << "Grad() " << ij << "=>" << fp << endl;})
					bool isInside = _I1->isInside(fp[0], fp[1])
							&& _I1->isInside(fr1[0], fr1[1])
							&& _I1->isInside(fr0[0], fr0[1])
							&& _I1->isInside(fp[0] + 1, fp[0])
							&& _I1->isInside(fp[0] - 1, fp[0]);
					if (isInside) {
						int I0ij = I0.bilinear(i, j);
						int I1xy = I1.bilinear(fp[0], fp[1]);
						int f = I0ij - I1xy;
						dtx += -f
								* (I1.bilinear(fp[0] + 1, fp[1])
										- I1.bilinear(fp[0] - 1, fp[1])) / 2.f;
						dty += -f
								* (I1.bilinear(fp[0], fp[1] + 1)
										- I1.bilinear(fp[0], fp[1] - 1)) / 2.f;
						dAngle += -f
								* (I1.bilinear(fr1[0], fr1[1])
										- I1.bilinear(fr0[0], fr0[1])) / 2.f;
						dScale += -f
								* (I1.bilinear(fs1[0], fs1[1])
										- I1.bilinear(fs0[0], fs0[1])) / 2.f;
						validCount++;
					}
				}
			}
		}
		validCount = 1;
		dtx /= validCount;
		dty /= validCount;
		dAngle /= validCount;
		dScale /= validCount;


		grad.tx = dtx;
		grad.ty = dty;
		grad.theta = 0; //dAngle;
	}

	void hessian(SquareMatrixR<float, 8>& sq) {
	}
};

