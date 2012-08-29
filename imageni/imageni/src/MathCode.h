/*
 * CMath.h
 *
 *  Created on: Jul 20, 2012
 *      Author: joohwile
 */

#ifndef CMATH_H_
#define CMATH_H_

#include <string>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include "Macros.h"

using namespace std;

namespace cmath {

const float EPS = 1e-5;
const float PI = 3.141592653589793238462643383;
const float PI_2 = 3.141592653589793238462643383 / 2.0;

float abs(float n);
int round(float n);
float sqrt(float s);
float cos(float x);
float sin(float e);

template<typename ElemType, int N>
class CVec {
public:
	ElemType _x[N];
	int length;
	CVec() :
			length(N) {
		zero();
	}
	CVec(ElemType x) :
			length(N) {
		zero();
		_x[0] = x;
	}
	CVec(ElemType x, ElemType y) :
			length(N) {
		zero();
		_x[0] = x;
		_x[1] = y;
	}
	CVec(ElemType x, ElemType y, ElemType z) :
			length(N) {
		zero();
		_x[0] = x;
		_x[1] = y;
		_x[2] = z;
	}
	CVec(ElemType x, ElemType y, ElemType z, ElemType w) :
			length(N) {
		zero();
		_x[0] = x;
		_x[1] = y;
		_x[2] = z;
		_x[3] = w;
	}
	int n() {
		return N;
	}
	void zero() {
		memset(_x, 0, sizeof(_x));
	}
	ElemType& operator[](const int j) {
		return _x[j];
	}
	ElemType operator[](const int j) const {
		return _x[j];
	}
	inline ElemType& x() {
		return _x[0];
	}
	inline ElemType& y() {
		return _x[1];
	}
	inline ElemType& z() {
		return _x[2];
	}
	inline ElemType x() const {
		return _x[0];
	}
	inline ElemType y() const {
		return _x[1];
	}
	inline ElemType z() const {
		return _x[2];
	}
	CVec<ElemType, N> operator+(const CVec<ElemType, N>& a) const {
		CVec<ElemType, N> o;
		for (int i = 0; i < N; i++) {
			o._x[i] = _x[i] + a._x[i];
		}
		return o;
	}
	CVec<ElemType, N> operator-(const CVec<ElemType, N>& a) const {
		CVec<ElemType, N> o;
		for (int i = 0; i < N; i++) {
			o._x[i] = _x[i] - a._x[i];
		}
		return o;
	}
	CVec<ElemType, N> operator*(const ElemType a) const {
		CVec<ElemType, N> o;
		for (int i = 0; i < N; i++) {
			o._x[i] = a * _x[i];
		}
		return o;
	}
	CVec<ElemType, N> operator*(const CVec<ElemType, N>& a) const {
		CVec<ElemType, N> o;
		for (int i = 0; i < N; i++) {
			o._x[i] = _x[i] * a._x[i];
		}
		return o;
	}
	CVec<ElemType, N> operator/(const ElemType b) const {
		CVec<ElemType, N> o;
		if (b == 0) {
			// treat as division by 1
			return (*this);
		}
		for (int i = 0; i < N; i++) {
			o._x[i] = _x[i] / b;
		}
		return o;
	}
	void negate() {
		forN(i)
		{
			_x[i] = -_x[i];
		}
	}
	void homogenize() {
		_x[0] /= _x[N - 1];
		_x[1] /= _x[N - 1];
	}
	void mult(const ElemType a, CVec<ElemType, N>& o) {
		forN(i)
		{
			o._x[i] = a * _x[i];
		}
	}
	void plus(const CVec<ElemType, N>& a, CVec<ElemType, N>& o) {
		forN(i)
		{
			o._x[i] = a._x[i] + _x[i];
		}
	}
	void linearCombination(const ElemType a, const ElemType b,
			const CVec<ElemType, N>& y, CVec<ElemType, N>& o) {
		forN(i)
		{
			o._x[i] = a * _x[i] + b * y._x[i];
		}
	}
	CVec<ElemType, N> round() const {
		CVec<ElemType, N> o;
		for (int i = 0; i < N; i++) {
			o._x[i] = cmath::round(_x[i]);
		}
		return o;
	}
	CVec<ElemType, N> pow2() const {
		CVec<ElemType, N> o;
		for (int i = 0; i < N; i++) {
			o._x[i] = _x[i] * _x[i];
		}
		return o;
	}
	ElemType sum() const {
		ElemType s = 0;
		for (int i = 0; i < N; i++) {
			s += _x[i];
		}
		return s;
	}
	ElemType squaredSum() const {
		ElemType s = 0;
		for (int i = 0; i < N; i++) {
			s += (_x[i] * _x[i]);
		}
		return s;
	}
	ElemType mag() const {
		ElemType s = squaredSum();
		return sqrt(s);
	}
	CVec<ElemType, N> normal() const {
		ElemType m = mag();
		CVec<ElemType, N> n;
		if (m == 0) {
			n.zero();
			return n;
		}
		for (int i = 0; i < N; i++) {
			n._x[i] = (_x[i] / m);
		}
		return n;
	}
	CVec<ElemType, N> neighbor(int dir, ElemType dist) const {
		CVec<ElemType, N> o;
		memcpy(o._x, _x, sizeof(_x));
		o[dir] += dist;
		return o;
	}
	CVec<ElemType, N>& set(ElemType a) {
		_x[0] = a;
		return (*this);
	}
	CVec<ElemType, N>& set(ElemType a, ElemType b) {
		_x[0] = a;
		_x[1] = b;
		return (*this);
	}
	CVec<ElemType, N>& point(ElemType a, ElemType b) {
		_x[0] = a;
		_x[1] = b;
		if (N > 2) {
			_x[2] = 1;
		}
		return (*this);
	}
	CVec<ElemType, N>& set(ElemType a, ElemType b, ElemType c) {
		_x[0] = a;
		_x[1] = b;
		_x[2] = c;
		return (*this);
	}
	CVec<ElemType, N>& point(ElemType a, ElemType b, ElemType c) {
		_x[0] = a;
		_x[1] = b;
		_x[2] = c;
		if (N > 3) {
			_x[3] = 1;
		}
		return (*this);
	}
	CVec<ElemType, N>& set(ElemType a, ElemType b, ElemType c, ElemType d) {
		_x[0] = a;
		_x[1] = b;
		_x[2] = c;
		_x[3] = d;
		return (*this);
	}
	CVec<ElemType, N>& setN(ElemType* a) {
		for (int i = 0; i < N; i++) {
			_x[i] = a[i];
		}
		return (*this);
	}
	CVec<ElemType, N>& setN(ElemType* a, int n) {
		for (int i = 0; i < n; i++) {
			_x[i] = a[i];
		}
		return (*this);
	}
};

typedef CVec<float, 2> Vec2;
typedef CVec<float, 3> Vec3;
typedef CVec<float, 4> Vec4;
typedef CVec<double, 4> Vec4d;

template<typename ElemType, int N>
inline CVec<ElemType, N> operator*(const float a, const CVec<ElemType, N>& v) {
	return (v * a);
}

template<typename ElemType, int N>
ostream& operator<<(ostream& o, CVec<ElemType, N> v) {
	//o << "(" << v[0];
	o << v[0];
	for (int i = 1; i < N; i++) {
		o << "," << v[i];
	}
	//o << ")";
	return o;
}

static const float _eps = 1e-6f;

template<typename T, int N>
class SquareMatrixR {
private:

public:
	T _d[N * N];
	T* _r[N];

	SquareMatrixR() {
		for (int i = 0; i < N; i++) {
			_r[i] = &_d[i * N];
		}
		identity();
	}
	inline void copyFrom(SquareMatrixR<T, N>& in) {
		for (int i = 0; i < size(); i++) {
			_d[i] = in._d[i];
		}
	}
	int size() const {
		return N * N;
	}
	inline SquareMatrixR<T, N>& set(int r, T a) {
		_r[r][0] = a;
		return (*this);
	}
	inline SquareMatrixR<T, N>& set(int r, T a, T b) {
		_r[r][0] = a;
		_r[r][1] = b;
		return (*this);
	}
	inline SquareMatrixR<T, N>& set(int r, T a, T b, T c) {
		_r[r][0] = a;
		_r[r][1] = b;
		_r[r][2] = c;
		return (*this);
	}
	inline SquareMatrixR<T, N>& set(int r, T a, T b, T c, T d) {
		_r[r][0] = a;
		_r[r][1] = b;
		_r[r][2] = c;
		_r[r][3] = d;
		return (*this);
	}
	inline void zero() {
		memset(_d, 0, sizeof(_d));
	}
	T& operator()(int i) {
		return _d[i];
	}
	T& operator()(int x, int y) {
		return _r[y][x];
	}
	T* operator[](int n) {
		return _r[n];
	}
	T det() {
		switch (N) {
		case 2:
			return _d[0] * _d[3] - _d[1] * _d[2];
		case 3:
			return _d[0] * (_d[4] * _d[8] - _d[5] * _d[7])
					+ _d[1] * (_d[5] * _d[6] - _d[3] * _d[8])
					+ _d[2] * (_d[3] * _d[7] - _d[4] * _d[6]);
		}
		return 0;
	}
	SquareMatrixR<T, N>& inverse(SquareMatrixR<T, N>& o, T* det) {
		if (N == 2) {
			o._d[0] = _d[3];
			o._d[1] = -_d[1];
			o._d[2] = -_d[2];
			o._d[3] = _d[0];
			*det = _d[0] * _d[3] - _d[1] * _d[2];
		} else if (N == 3) {
			o._d[0] = _d[4] * _d[8] - _d[5] * _d[7];
			o._d[1] = _d[7] * _d[2] - _d[1] * _d[8];
			o._d[2] = _d[1] * _d[5] - _d[2] * _d[4];
			o._d[3] = _d[5] * _d[6] - _d[3] * _d[8];
			o._d[4] = _d[0] * _d[8] - _d[2] * _d[6];
			o._d[5] = _d[2] * _d[3] - _d[5] * _d[0];
			o._d[6] = _d[3] * _d[7] - _d[4] * _d[6];
			o._d[7] = _d[1] * _d[6] - _d[7] * _d[0];
			o._d[8] = _d[0] * _d[4] - _d[1] * _d[3];
			*det = _d[0] * o._d[0] + _d[1] * o._d[3] + _d[2] * o._d[6];
		} else if (N == 4) {
			o._d[0] = _d[5] * _d[10] * _d[15] + _d[6] * _d[11] * _d[13]
					+ _d[7] * _d[9] * _d[14] - _d[5] * _d[11] * _d[14]
					- _d[6] * _d[9] * _d[15] - _d[7] * _d[10] * _d[13];
			o._d[1] = _d[1] * _d[11] * _d[14] + _d[2] * _d[9] * _d[15]
					+ _d[3] * _d[10] * _d[13] - _d[1] * _d[10] * _d[15]
					- _d[2] * _d[11] * _d[13] - _d[3] * _d[9] * _d[14];
			o._d[2] = _d[1] * _d[6] * _d[15] + _d[2] * _d[7] * _d[13]
					+ _d[3] * _d[5] * _d[14] - _d[1] * _d[7] * _d[14]
					- _d[2] * _d[5] * _d[15] - _d[3] * _d[6] * _d[13];
			o._d[3] = _d[1] * _d[7] * _d[10] + _d[2] * _d[5] * _d[11]
					+ _d[3] * _d[6] * _d[9] - _d[1] * _d[6] * _d[11]
					- _d[2] * _d[7] * _d[9] - _d[3] * _d[5] * _d[10];
			o._d[4] = _d[4] * _d[11] * _d[14] + _d[6] * _d[8] * _d[15]
					+ _d[7] * _d[10] * _d[12] - _d[4] * _d[10] * _d[15]
					- _d[6] * _d[11] * _d[12] - _d[7] * _d[8] * _d[14];
			o._d[5] = _d[0] * _d[10] * _d[15] + _d[2] * _d[11] * _d[12]
					+ _d[3] * _d[8] * _d[14] - _d[0] * _d[11] * _d[14]
					- _d[2] * _d[8] * _d[15] - _d[3] * _d[10] * _d[12];
			o._d[6] = _d[0] * _d[7] * _d[14] + _d[2] * _d[4] * _d[15]
					+ _d[3] * _d[6] * _d[12] - _d[0] * _d[6] * _d[15]
					- _d[2] * _d[7] * _d[12] - _d[3] * _d[4] * _d[14];
			o._d[7] = _d[0] * _d[6] * _d[11] + _d[2] * _d[7] * _d[8]
					+ _d[3] * _d[4] * _d[10] - _d[0] * _d[7] * _d[10]
					- _d[2] * _d[4] * _d[11] - _d[3] * _d[6] * _d[8];
			o._d[8] = _d[4] * _d[9] * _d[15] + _d[5] * _d[11] * _d[12]
					+ _d[7] * _d[8] * _d[13] - _d[4] * _d[11] * _d[13]
					- _d[5] * _d[8] * _d[15] - _d[7] * _d[9] * _d[12];
			o._d[9] = _d[0] * _d[11] * _d[13] + _d[1] * _d[8] * _d[15]
					+ _d[3] * _d[9] * _d[12] - _d[0] * _d[9] * _d[15]
					- _d[1] * _d[11] * _d[12] - _d[3] * _d[8] * _d[13];
			o._d[10] = _d[0] * _d[5] * _d[15] + _d[1] * _d[7] * _d[12]
					+ _d[3] * _d[4] * _d[13] - _d[0] * _d[7] * _d[13]
					- _d[1] * _d[4] * _d[15] - _d[3] * _d[5] * _d[12];
			o._d[11] = _d[0] * _d[7] * _d[9] + _d[1] * _d[4] * _d[11]
					+ _d[3] * _d[5] * _d[8] - _d[0] * _d[5] * _d[11]
					- _d[1] * _d[7] * _d[8] - _d[3] * _d[4] * _d[9];
			o._d[12] = _d[4] * _d[10] * _d[13] + _d[5] * _d[8] * _d[14]
					+ _d[6] * _d[9] * _d[12] - _d[4] * _d[9] * _d[14]
					- _d[5] * _d[10] * _d[12] - _d[6] * _d[8] * _d[13];
			o._d[13] = _d[0] * _d[9] * _d[14] + _d[1] * _d[10] * _d[12]
					+ _d[2] * _d[8] * _d[13] - _d[0] * _d[10] * _d[13]
					- _d[1] * _d[8] * _d[14] - _d[2] * _d[9] * _d[12];
			o._d[14] = _d[0] * _d[6] * _d[13] + _d[1] * _d[4] * _d[14]
					+ _d[2] * _d[5] * _d[12] - _d[0] * _d[5] * _d[14]
					- _d[1] * _d[6] * _d[12] - _d[2] * _d[4] * _d[13];
			o._d[15] = _d[0] * _d[5] * _d[10] + _d[1] * _d[6] * _d[8]
					+ _d[2] * _d[4] * _d[9] - _d[0] * _d[6] * _d[9]
					- _d[1] * _d[4] * _d[10] - _d[2] * _d[5] * _d[8];

			*det = _d[0] * o._d[0] + _d[1] * o._d[4] + _d[2] * o._d[8]
					+ _d[3] * o._d[12];
			/*
			 T* mat = _d;
			 T* dst = o._d;
			 T tmp[12]; // temp array for pairs
			 T src[16]; // array of transpose source matrix
			 // transpose matrix
			 for (int i = 0; i < 4; i++) {
			 src[i] = mat[i * 4];
			 src[i + 4] = mat[i * 4 + 1];
			 src[i + 8] = mat[i * 4 + 2];
			 src[i + 12] = mat[i * 4 + 3];
			 }
			 // calculate pairs for first 8 elements (cofactors)
			 tmp[0] = src[10] * src[15];
			 tmp[1] = src[11] * src[14];
			 tmp[2] = src[9] * src[15];
			 tmp[3] = src[11] * src[13];
			 tmp[4] = src[9] * src[14];
			 tmp[5] = src[10] * src[13];
			 tmp[6] = src[8] * src[15];
			 tmp[7] = src[11] * src[12];
			 tmp[8] = src[8] * src[14];
			 tmp[9] = src[10] * src[12];
			 tmp[10] = src[8] * src[13];
			 tmp[11] = src[9] * src[12];
			 // calculate first 8 elements (cofactors)
			 dst[0] = tmp[0] * src[5] + tmp[3] * src[6] + tmp[4] * src[7];
			 dst[0] -= tmp[1] * src[5] + tmp[2] * src[6] + tmp[5] * src[7];
			 dst[1] = tmp[1] * src[4] + tmp[6] * src[6] + tmp[9] * src[7];
			 dst[1] -= tmp[0] * src[4] + tmp[7] * src[6] + tmp[8] * src[7];
			 dst[2] = tmp[2] * src[4] + tmp[7] * src[5] + tmp[10] * src[7];
			 dst[2] -= tmp[3] * src[4] + tmp[6] * src[5] + tmp[11] * src[7];
			 dst[3] = tmp[5] * src[4] + tmp[8] * src[5] + tmp[11] * src[6];
			 dst[3] -= tmp[4] * src[4] + tmp[9] * src[5] + tmp[10] * src[6];
			 dst[4] = tmp[1] * src[1] + tmp[2] * src[2] + tmp[5] * src[3];
			 dst[4] -= tmp[0] * src[1] + tmp[3] * src[2] + tmp[4] * src[3];
			 dst[5] = tmp[0] * src[0] + tmp[7] * src[2] + tmp[8] * src[3];
			 dst[5] -= tmp[1] * src[0] + tmp[6] * src[2] + tmp[9] * src[3];
			 dst[6] = tmp[3] * src[0] + tmp[6] * src[1] + tmp[11] * src[3];
			 dst[6] -= tmp[2] * src[0] + tmp[7] * src[1] + tmp[10] * src[3];
			 dst[7] = tmp[4] * src[0] + tmp[9] * src[1] + tmp[10] * src[2];
			 dst[7] -= tmp[5] * src[0] + tmp[8] * src[1] + tmp[11] * src[2];
			 // calculate pairs for second 8 elements (cofactors)
			 tmp[0] = src[2] * src[7];
			 tmp[1] = src[3] * src[6];
			 tmp[2] = src[1] * src[7];
			 tmp[3] = src[3] * src[5];
			 tmp[4] = src[1] * src[6];
			 tmp[5] = src[2] * src[5];
			 tmp[6] = src[0] * src[7];
			 tmp[7] = src[3] * src[4];
			 tmp[8] = src[0] * src[6];
			 tmp[9] = src[2] * src[4];
			 tmp[10] = src[0] * src[5];
			 tmp[11] = src[1] * src[4];
			 // calculate second 8 elements (cofactors)
			 dst[8] = tmp[0] * src[13] + tmp[3] * src[14] + tmp[4] * src[15];
			 dst[8] -= tmp[1] * src[13] + tmp[2] * src[14] + tmp[5] * src[15];
			 dst[9] = tmp[1] * src[12] + tmp[6] * src[14] + tmp[9] * src[15];
			 dst[9] -= tmp[0] * src[12] + tmp[7] * src[14] + tmp[8] * src[15];
			 dst[10] = tmp[2] * src[12] + tmp[7] * src[13] + tmp[10] * src[15];
			 dst[10] -= tmp[3] * src[12] + tmp[6] * src[13] + tmp[11] * src[15];
			 dst[11] = tmp[5] * src[12] + tmp[8] * src[13] + tmp[11] * src[14];
			 dst[11] -= tmp[4] * src[12] + tmp[9] * src[13] + tmp[10] * src[14];
			 dst[12] = tmp[2] * src[10] + tmp[5] * src[11] + tmp[1] * src[9];
			 dst[12] -= tmp[4] * src[11] + tmp[0] * src[9] + tmp[3] * src[10];
			 dst[13] = tmp[8] * src[11] + tmp[0] * src[8] + tmp[7] * src[10];
			 dst[13] -= tmp[6] * src[10] + tmp[9] * src[11] + tmp[1] * src[8];
			 dst[14] = tmp[6] * src[9] + tmp[11] * src[11] + tmp[3] * src[8];
			 dst[14] -= tmp[10] * src[11] + tmp[2] * src[8] + tmp[7] * src[9];
			 dst[15] = tmp[10] * src[10] + tmp[4] * src[8] + tmp[9] * src[9];
			 dst[15] -= tmp[8] * src[9] + tmp[11] * src[10] + tmp[5] * src[8];
			 // calculate determinant
			 *det = src[0] * dst[0] + src[1] * dst[1] + src[2] * dst[2] + src[3] * dst[3];
			 // calculate matrix inverse
			 double d = 1 / *det;
			 for (int j = 0; j < 16; j++) {
			 dst[j] *= d;
			 }
			 */
		}
		return o;
	}
	SquareMatrixR<T, N>& inverse(SquareMatrixR<T, N>& o) {
		T det;
		inverse(o, &det);
		o.divide(det);
		return (o);
	}

	void suppressEpsilion() {
		for (int i = 0; i < N * N; i++) {
			_d[i] = (_d[i] < _eps && _d[i] > -_eps) ? 0.0 : _d[i];
		}
	}
	void identity() {
		zero();
		for (int i = 0; i < N; i++) {
			_r[i][i] = 1;
		}
	}
	SquareMatrixR<T, N>& transpose() {
		for (int j = 0; j < N; j++) {
			for (int i = 0; i < N; i++) {
				if (i < j) {
					T t = _r[i][j];
					_r[j][i] = _r[i][j];
					_r[j][i] = t;
				}
			}
		}
		return (*this);
	}
	SquareMatrixR<T, N>& transpose(SquareMatrixR<T, N>& c) {
		for (int j = 0; j < N; j++) {
			for (int i = 0; i < N; i++) {
				c._r[j][i] = _r[i][j];
			}
		}
		return c;
	}

	SquareMatrixR<T, N>& image2standard(int w, int h) {
		if (N == 3) {
			//_d[0] = -_d[0];
			//_d[2] = _d[2] + w;
			_d[4] = -_d[4];
			_d[5] = _d[5] + h;
		}
		return (*this);
	}
	SquareMatrixR<T, N>& standard2image(int w, int h) {
		if (N == 3) {
			_d[4] = -_d[4];
			_d[5] = _d[5] - h;
		}
		return (*this);
	}
	/*
	 CSquareMat<N>& rotation2D(T deg, T cx, T cy) {
	 T rad = deg / 180 * PI;
	 T cr = cos(rad), sr = sin(rad);
	 zero();
	 _d[0][0] = cr;
	 _d[0][1] = -sr;
	 _d[1][0] = sr;
	 _d[1][1] = cr;
	 _d[2][2] = 1;
	 return (*this);
	 }
	 *
	 */

	template<int M>
	SquareMatrixR<T, N>& scale(CVec<T, M>& scale) {
		if (M < N) {
			for (int i = 0; i < M; i++) {
				_r[i][i] = scale[i] * _r[i][i];
			}
		}
		return (*this);
	}

	SquareMatrixR<T, N>& rotateAxisAround(T deg, T cx, T cy) {
		if (deg < _eps && deg > _eps) {
			return (*this);
		}
		T rad = deg / 180 * PI;
		T cr = cos(rad), sr = sin(rad);
		if (N == 3) {
			T a = _d[0];
			T b = _d[1];
			T c = _d[2] - cx;
			T d = _d[3];
			T e = _d[4];
			T f = _d[5] - cy;
			_d[0] = a * cr + d * sr;
			_d[1] = b * cr + e * sr;
			_d[2] = c * cr + f * sr + cx;
			_d[3] = a * -sr + d * cr;
			_d[4] = b * -sr + e * cr;
			_d[5] = c * -sr + f * cr + cy;
		} else if (N == 2) {
			T a = _d[0];
			T b = _d[1];
			T c = _d[2];
			T d = _d[3];
			_d[0] = cr * a + sr * c;
			_d[1] = cr * b + sr * d;
			_d[2] = -sr * a + cr * c;
			_d[3] = -sr * b + cr * d;
		}
		return (*this);
	}
	SquareMatrixR<T, N>& scale(Vec3& v) {
		_r[0][0] *= v[0];
		_r[1][1] *= v[1];
		_r[2][2] *= v[2];
		return (*this);
	}
	SquareMatrixR<T, N>& translateOrigin(Vec2& v) {
		_r[0][N - 1] = _r[0][N - 1] - v[0];
		_r[1][N - 1] = _r[1][N - 1] - v[1];
		return (*this);
	}
	SquareMatrixR<T, N>& translateOrigin(double tx, double ty) {
		_r[0][N - 1] = _r[0][N - 1] - tx;
		_r[1][N - 1] = _r[1][N - 1] - ty;
		return (*this);
	}

	SquareMatrixR<T, N>& translateOrigin(Vec3& v) {
		_r[0][N - 1] = _r[0][N - 1] - v[0];
		_r[1][N - 1] = _r[1][N - 1] - v[1];
		_r[2][N - 1] = _r[2][N - 1] - v[2];
		return (*this);
	}

	SquareMatrixR<T, N>& translateOriginN(CVec<T, N>& v) {
		for (int i = 0; i < N - 1; i++) {
			_r[i][N - 1] = _r[i][N - 1] - v[i];
		}
		return (*this);
	}

	SquareMatrixR<T, N>& plus(const T c) {
		for (int i = 0; i < N * N; i++) {
			_d[i] = _d[i] + c;
		}
		return (*this);
	}
	SquareMatrixR<T, N>& minus(T c) {
		for (int i = 0; i < N * N; i++) {
			_d[i] = _d[i] - c;
		}
		return (*this);
	}
	SquareMatrixR<T, N>& mult(T c) {
		for (int i = 0; i < N * N; i++) {
			_d[i] = _d[i] * c;
		}
		return (*this);
	}
	SquareMatrixR<T, N>& divide(T c) {
		for (int i = 0; i < N * N; i++) {
			_d[i] = _d[i] / c;
		}
		return (*this);
	}
	SquareMatrixR<T, N>& plus(SquareMatrixR<T, N>&b, SquareMatrixR<T, N>&c) {
		for (int i = 0; i < N * N; i++) {
			c._d[i] = _d[i] + b._d[i];
		}
		return c;
	}
	SquareMatrixR<T, N>& minus(SquareMatrixR<T, N>&b, SquareMatrixR<T, N>&c) {
		for (int i = 0; i < N * N; i++) {
			c._d[i] = _d[i] + b._d[i];
		}
		return c;
	}
	SquareMatrixR<T, N>& dot(SquareMatrixR<T, N>&b, SquareMatrixR<T, N>&c) {
		for (int i = 0; i < N * N; i++) {
			c._d[i] = _d[i] * b._d[i];
		}
		return c;
	}
	SquareMatrixR<T, N>& divide(SquareMatrixR<T, N>&b, SquareMatrixR<T, N>&c) {
		for (int i = 0; i < N * N; i++) {
			c._d[i] = _d[i] / b._d[i];
		}
		return c;
	}
	T* mult(T* b, T* c) {
		for (int i = 0; i < N; i++) {
			T e = 0;
			for (int l = 0; l < N; l++) {
				e += _r[i][l] * b[l];
			}
			c[i] = e;
		}
		return c;
	}
	template<typename U>
	CVec<T, N>& mult(U* b, CVec<T, N> &c) {
		for (int i = 0; i < N; i++) {
			U e = 0;
			for (int l = 0; l < N; l++) {
				e += _r[i][l] * b[l];
			}
			c[i] = e;
		}
		return c;
	}

	template<typename U>
	CVec<T, N>& mult(CVec<T, N>& b, CVec<T, N>& c) {
		for (int i = 0; i < N; i++) {
			U e = 0;
			for (int l = 0; l < N; l++) {
				e += _r[i][l] * b[l];
			}
			c[i] = e;
		}
		return c;
	}
	SquareMatrixR<T, N>& mult(SquareMatrixR<T, N>& b, SquareMatrixR<T, N>& c) {
		for (int j = 0; j < N; j++) {
			for (int i = 0; i < N; i++) {
				T e = 0;
				for (int l = 0; l < N; l++) {
					e += (_r[j][l]) * b[l][i];
				}
				c[j][i] = e;
			}
		}
		return c;
	}
	friend ostream& operator<<(ostream& os, SquareMatrixR<T, N> &m) {
		os << "[";
		for (int j = 0; j < N; j++) {
			if (j > 0) {
				os << "; ";
			}
			os << m[j][0];
			for (int i = 1; i < N; i++) {
				os << ",";
				os << m[j][i];
			}
		}
		os << "]";
		os.flush();
		return os;
	}
};

template<typename T>
class MatrixR {
public:
	int nRows;
	int nCols;
	int nElems;
	T* _M;
	T** _R;
public:
	MatrixR(int m, int n) {
		nRows = m;
		nCols = n;
		nElems = m * n;
		_M = new T[m * n];
		_R = new T*[m];
		for (int i = 0; i < m; i++) {
			_R[i] = &_M[i * n];
		}
		zero();
	}
	virtual ~MatrixR() {
		delete[] _R;
		delete[] _M;
	}
	void zero() {
		for (int i = 0; i < nElems; i++) {
			_M[i] = 0;
		}
	}
	T& operator()(int r, int c) {
		return _R[r][c];
	}
	T* operator[](int n) {
		return _R[n];
	}
	void set(T* row, int r) {
		for (int i = 0; i < nCols; i++) {
			_R[r][i] = row[i];
		}
	}
	void set(T a, T b, T c, T d, int r) {
		_R[r][0] = a;
		_R[r][1] = b;
		_R[r][2] = c;
		_R[r][3] = d;
	}
	void transposeInto(MatrixR<T>& o) {
		for (int i = 0; i < o.nRows; i++) {
			for (int j = 0; j < o.nCols; j++) {
				o[i][j] = _R[j][i];
			}
		}
	}
};

typedef MatrixR<float> DynMatF;
typedef MatrixR<double> DynMatD;

template<typename T>
std::ostream& operator<<(std::ostream& os, const MatrixR<T>& mat) {
	os << "Matrix: " << mat.nRows << " x " << mat.nCols << endl;
	for (int i = 0; i < __cmath_min(mat.nRows, 3); i++) {
		os << mat._R[i][0];
		for (int j = 1; j < __cmath_min(mat.nCols, 10); j++) {
			os << ", " << mat._R[i][j];
		}
		os << endl;
	}
	return os;
}

template<typename T>
inline T dot(MatrixR<T> &a, int row, MatrixR<T> &b, int col);

template<typename T>
void mult(MatrixR<T>& m, MatrixR<T>& n, MatrixR<T>& o);

typedef SquareMatrixR<float, 2> Mat2;
typedef SquareMatrixR<float, 3> Mat3;
typedef SquareMatrixR<float, 4> Mat4;
typedef SquareMatrixR<double, 4> Mat4d;

Vec4 vecmin(Vec4* a, int n);
Vec3 vecmin(Vec3* a, int n);
Vec2 vecmin(Vec2* a, int n);
Vec4 vecmax(Vec4* a, int n);
Vec3 vecmax(Vec3* a, int n);
Vec2 vecmax(Vec2* a, int n);

template<typename T, int N>
CVec<T, N>& minN(CVec<T, N>* a, int n, CVec<T, N>& o) {
	o = a[0];
	for (int i = 1; i < n; i++) {
		for (int j = 0; j < N; j++) {
			o._x[j] = o._x[j] < a[i]._x[j] ? a[i]._x[j] : o._x[j];
		}
	}
	return o;
}

inline Vec2& mult(Mat3& a, Vec2& b, Vec2& c) {
	c[0] = a._d[0] * b[0] + a._d[1] * b[1] + a._d[2];
	c[1] = a._d[3] * b[0] + a._d[4] * b[1] + a._d[5];
	return c;
}

inline Vec3& mult(Mat3& a, Vec3& b, Vec3& c) {
	c._x[0] = a._d[0] * b._x[0] + a._d[1] * b._x[1] + a._d[2] * b._x[2];
	c._x[1] = a._d[3] * b._x[0] + a._d[4] * b._x[1] + a._d[5] * b._x[2];
	c._x[2] = a._d[6] * b._x[0] + a._d[7] * b._x[1] + a._d[8] * b._x[2];
	return c;
}

inline void scaleZ(Vec3* a, int n, Vec3 *b) {
	for (int i = 0; i < n; i++) {
		if (abs(a[i]._x[2]) > cmath::EPS) {
			b[i]._x[0] = a[i]._x[0] / a[i]._x[2];
			b[i]._x[1] = a[i]._x[1] / a[i]._x[2];
			b[i]._x[2] = a[i]._x[2] / a[i]._x[2];
		}
	}
}

inline void mult(Mat3& a, Vec2* b, Vec2* c, int n) {
	for (int i = 0; i < n; i++) {
		mult(a, b[i], c[i]);
	}
}

inline void mult(Mat3& a, Vec3* b, Vec3* c, int n) {
	for (int i = 0; i < n; i++) {
		mult(a, b[i], c[i]);
	}
}

inline Vec4& mult(Mat4& a, Vec4& b, Vec4& c) {
	c[0] = a._d[0] * b[0] + a._d[1] * b[1] + a._d[2] * b[2] + a._d[3] * b[3];
	c[1] = a._d[4] * b[0] + a._d[5] * b[1] + a._d[6] * b[2] + a._d[7] * b[3];
	c[2] = a._d[8] * b[0] + a._d[9] * b[1] + a._d[10] * b[2] + a._d[11] * b[3];
	c[3] = a._d[12] * b[0] + a._d[13] * b[1] + a._d[14] * b[2]
			+ a._d[15] * b[3];
	return c;
}

inline Vec4d& mult(Mat4d& a, Vec4d& b, Vec4d& c) {
	c[0] = a._d[0] * b[0] + a._d[1] * b[1] + a._d[2] * b[2] + a._d[3] * b[3];
	c[1] = a._d[4] * b[0] + a._d[5] * b[1] + a._d[6] * b[2] + a._d[7] * b[3];
	c[2] = a._d[8] * b[0] + a._d[9] * b[1] + a._d[10] * b[2] + a._d[11] * b[3];
	c[3] = a._d[12] * b[0] + a._d[13] * b[1] + a._d[14] * b[2]
			+ a._d[15] * b[3];
	return c;
}

inline void mult(Mat4& a, Vec4* b, Vec4* c, int n) {
	for (int i = 0; i < n; i++) {
		mult(a, b[i], c[i]);
	}
}

template<typename T>
inline SquareMatrixR<T, 2>& mult(SquareMatrixR<T, 2>& a, SquareMatrixR<T, 2>& b,
		SquareMatrixR<T, 2>& c) {
	T e;
	e = 0;
	e += a[0][0] * b[0][0];
	e += a[0][1] * b[1][0];
	c[0][0] = e;
	e = 0;
	e += a[0][0] * b[0][1];
	e += a[0][1] * b[1][1];
	c[0][1] = e;
	e = 0;
	e += a[1][0] * b[0][0];
	e += a[1][1] * b[1][0];
	c[1][0] = e;
	e = 0;
	e += a[1][0] * b[0][1];
	e += a[1][1] * b[1][1];
	c[1][1] = e;
	return c;
}

template<typename T>
inline SquareMatrixR<T, 3> mult(SquareMatrixR<T, 3>& a, SquareMatrixR<T, 3>& b,
		SquareMatrixR<T, 3>& c) {
	T e;
	e = 0;
	e += a[0][0] * b[0][0];
	e += a[0][1] * b[1][0];
	e += a[0][2] * b[2][0];
	c[0][0] = e;
	e = 0;
	e += a[0][0] * b[0][1];
	e += a[0][1] * b[1][1];
	e += a[0][2] * b[2][1];
	c[0][1] = e;
	e = 0;
	e += a[0][0] * b[0][2];
	e += a[0][1] * b[1][2];
	e += a[0][2] * b[2][2];
	c[0][2] = e;
	e = 0;
	e += a[1][0] * b[0][0];
	e += a[1][1] * b[1][0];
	e += a[1][2] * b[2][0];
	c[1][0] = e;
	e = 0;
	e += a[1][0] * b[0][1];
	e += a[1][1] * b[1][1];
	e += a[1][2] * b[2][1];
	c[1][1] = e;
	e = 0;
	e += a[1][0] * b[0][2];
	e += a[1][1] * b[1][2];
	e += a[1][2] * b[2][2];
	c[1][2] = e;
	e = 0;
	e += a[2][0] * b[0][0];
	e += a[2][1] * b[1][0];
	e += a[2][2] * b[2][0];
	c[2][0] = e;
	e = 0;
	e += a[2][0] * b[0][1];
	e += a[2][1] * b[1][1];
	e += a[2][2] * b[2][1];
	c[2][1] = e;
	e = 0;
	e += a[2][0] * b[0][2];
	e += a[2][1] * b[1][2];
	e += a[2][2] * b[2][2];
	c[2][2] = e;
	return c;
}

template<typename T>
inline SquareMatrixR<T, 4>& mult(SquareMatrixR<T, 4>& a, SquareMatrixR<T, 4>& b,
		SquareMatrixR<T, 4>& c) {
	//Mat4 c;
	T e = 0;
	e = 0;
	e += a[0][0] * b[0][0];
	e += a[0][1] * b[1][0];
	e += a[0][2] * b[2][0];
	e += a[0][3] * b[3][0];
	c[0][0] = e;
	e = 0;
	e += a[0][0] * b[0][1];
	e += a[0][1] * b[1][1];
	e += a[0][2] * b[2][1];
	e += a[0][3] * b[3][1];
	c[0][1] = e;
	e = 0;
	e += a[0][0] * b[0][2];
	e += a[0][1] * b[1][2];
	e += a[0][2] * b[2][2];
	e += a[0][3] * b[3][2];
	c[0][2] = e;
	e = 0;
	e += a[0][0] * b[0][3];
	e += a[0][1] * b[1][3];
	e += a[0][2] * b[2][3];
	e += a[0][3] * b[3][3];
	c[0][3] = e;
	e = 0;
	e += a[1][0] * b[0][0];
	e += a[1][1] * b[1][0];
	e += a[1][2] * b[2][0];
	e += a[1][3] * b[3][0];
	c[1][0] = e;
	e = 0;
	e += a[1][0] * b[0][1];
	e += a[1][1] * b[1][1];
	e += a[1][2] * b[2][1];
	e += a[1][3] * b[3][1];
	c[1][1] = e;
	e = 0;
	e += a[1][0] * b[0][2];
	e += a[1][1] * b[1][2];
	e += a[1][2] * b[2][2];
	e += a[1][3] * b[3][2];
	c[1][2] = e;
	e = 0;
	e += a[1][0] * b[0][3];
	e += a[1][1] * b[1][3];
	e += a[1][2] * b[2][3];
	e += a[1][3] * b[3][3];
	c[1][3] = e;
	e = 0;
	e += a[2][0] * b[0][0];
	e += a[2][1] * b[1][0];
	e += a[2][2] * b[2][0];
	e += a[2][3] * b[3][0];
	c[2][0] = e;
	e = 0;
	e += a[2][0] * b[0][1];
	e += a[2][1] * b[1][1];
	e += a[2][2] * b[2][1];
	e += a[2][3] * b[3][1];
	c[2][1] = e;
	e = 0;
	e += a[2][0] * b[0][2];
	e += a[2][1] * b[1][2];
	e += a[2][2] * b[2][2];
	e += a[2][3] * b[3][2];
	c[2][2] = e;
	e = 0;
	e += a[2][0] * b[0][3];
	e += a[2][1] * b[1][3];
	e += a[2][2] * b[2][3];
	e += a[2][3] * b[3][3];
	c[2][3] = e;
	e = 0;
	e += a[3][0] * b[0][0];
	e += a[3][1] * b[1][0];
	e += a[3][2] * b[2][0];
	e += a[3][3] * b[3][0];
	c[3][0] = e;
	e = 0;
	e += a[3][0] * b[0][1];
	e += a[3][1] * b[1][1];
	e += a[3][2] * b[2][1];
	e += a[3][3] * b[3][1];
	c[3][1] = e;
	e = 0;
	e += a[3][0] * b[0][2];
	e += a[3][1] * b[1][2];
	e += a[3][2] * b[2][2];
	e += a[3][3] * b[3][2];
	c[3][2] = e;
	e = 0;
	e += a[3][0] * b[0][3];
	e += a[3][1] * b[1][3];
	e += a[3][2] * b[2][3];
	e += a[3][3] * b[3][3];
	c[3][3] = e;
	return c;
}

}
#endif /* CMATH_H_ */
