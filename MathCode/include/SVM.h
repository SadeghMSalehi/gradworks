//
//  SVM.h
//  MathCode
//
//  Created by Joohwi Lee on 7/30/12.
//
//

#ifndef MathCode_SVM_h
#define MathCode_SVM_h

#include "MatrixCode.h"
#include "OptiCode.h"
#include <iostream>

using namespace std;

namespace MathCode {

template<int N>
class SVMCostFunc: public CFunc<DVec<float, N>, float, DVec<float, N>, int> {
public:
	InputType* _x;

	virtual ~SVMCostFunc() {
	}
	virtual void position(InputType& x) {
		_x = &x;
	}
	virtual void value(ValueOutputType& value) {

	}
	virtual void gradient(GradientOutputType& grad) {

	}
	virtual void hessian(HessianOutputType& hessian) {

	}
};

template<int N>
class SVM {
public:
	typedef DVec<float, N> InputType;
	int _nData;
	InputType _w;
	float _bias;
	InputType *_x;
	float* _alpha;
	float* _yxij;
	int* _y;

	SVM(InputType* data, int* y, int n) {
		_nData = n;
		_x = new InputType[_nData];
		_alpha = new float[_nData];
		_yxij = new float[_nData * (_nData + 1) / 2];
		_y = new int[_nData];
		for (int i = 0; i < n; i++) {
			_x[i].copyFrom(data[n]._V);
			_alpha[i] = 0;
			_y = y[i];
		}
		for (int i = 0; i < n; i++) {
			for (int j = i; j < n; j++) {
				setYXij(i, j, _y[i] * _y[j] * _x[i].dot(_x[j]));
			}
		}
		_w.zero();
		_bias = 0;
	}
	virtual ~SVM() {
		delete[] _x;
		delete[] _alpha;
		delete[] _y;
		delete[] _yxij;
	}
	inline void setYXij(int i, int j, float v) {
		_yxij[((j * (j + 1)) >> 1) + i] = v;
	}
	inline float getYXij(int i, int j) {
		return _yxij[((j * (j + 1)) >> 1) + i];
	}
	void printYXij(ostream& os) {
		for (int j = 0; j < _nData; j++) {
			for (int i = 0; i < _nData; i++) {
				int ii = __mathcode_min(i,j);
				int jj = __mathcode_max(i,j);
				cout << _yxij[(jj * (jj + 1)) >> 1 + ii] << endl;
			}
		}
	}
	void setSupport(InputType& w) {
		_w.copyFrom(w);
	}
	void setBias(int b) {
		_bias = b;
	}
	float getMargin() {
		float w = _w.mag();
		if (w == 0) {
			return 0;
		} else {
			return 1 / w;
		}
	}
	bool evaluate() {
		int positiveCount = 0;
		for (int i = 0; i < _nData; i++) {
			int fx = computeFx(i);
			if (fx >= 1) {
				positiveCount++;
			}
		}
		return positiveCount == N;
	}
	float computeFx(int i) {
		return _y[i] * (_x[i].dot(_w) - _bias);
	}
	float computeCost() {
		return InputType::NormSquare(_w._V);
	}
	float computeLagrangianCost(float* alpha = _alpha) {
		float cost = 0;
		for (int i = 0; i < _nData; i++) {
			cost += alpha[i];
			for (int j = 0; j < _nData; j++) {
				cost -= .5f * (alpha[i] * alpha[j] * _y[i] * _y[j] * _w[i].dot(_w[j]));
			}
		}
		return cost;
	}

	void computeGradientOfLagrangianCost(float* gradient, float* alpha = _alpha) {
		for (int i = 0; i < _nData; i++) {
			gradient[i] = 1;
			for (int j = 0; j < i; j++) {
				gradient[i] -= 0.5f * alpha[j] * getYXij(j, i);
			}
			gradient[i] -= getYXij(i, i);
			for (int j = i + 1; j < _nData; j++) {
				gradient[i] -= 0.5f * alpha[j] * getYXij(i, j);
			}
		}
	}
};

}
#endif
