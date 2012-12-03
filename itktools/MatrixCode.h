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

#define DEBUG
#define _DV_ DVec<T,N>
#define forN(i) for(int i = 0; i < N; i++)
#define printVar(v) std::cout << #v ":" << v << std::endl
#define printColumn(M,j) std::cout << #M "[" << j << "] = [ "; forN(_i_) std::cout << M[_i_][j] << " "; std::cout << "]" << std::endl
#define printRow(M,j) std::cout << #M "[" << j << "] = [ "; forN(_i_) std::cout << M[j][_i_] << " "; std::cout << "]" << std::endl
#define __mathcode_min(x,y) (x<y?x:y)
#define __mathcode_max(x,y) (x>y?x:y)
#define MAX_T 1e6
#define EPS_T 1e-5
#ifdef DEBUG
#define _DBG_(stmt) stmt
#endif

using namespace std;

namespace MathCode {

const float PI = 3.141592653589793238462643383;
const float PI_2 = 3.141592653589793238462643383 / 2.0;

float abs(float n);
int round(float n);
float sqrt(float s);
float cos(float x);
float sin(float e);

template<typename T, int N>
class DVec {
public:
	T _V[N];
	int length;

	DVec() {
		length = N;
		zero();
	}
	DVec(T* v) {
		length = N;
		forN(i)
		{
			_V[i] = v[i];
		}
	}
	DVec(DVec<T, N>* v) :
			length(N) {
		copyFrom(v->_V);
	}
	DVec(DVec<T, N>& v) :
			length(N) {
		copyFrom(v._V);
	}
	DVec(T a) :
			length(N) {
		_V[0] = a;
	}
	DVec(T a, T b) :
			length(N) {
		_V[0] = a;
		_V[1] = b;
	}
	DVec(T a, T b, T c) :
			length(N) {
		_V[0] = a;
		_V[1] = b;
		_V[2] = c;
	}
	DVec(T a, T b, T c, T d) :
			length(N) {
		_V[0] = a;
		_V[1] = b;
		_V[2] = c;
		_V[3] = d;
	}
	virtual ~DVec() {
	}
	static DVec<T, N>& New() {
		return *(new DVec());
	}
	static T NormSquare(T* v) {
		T n = 0;
		forN(i)
		{
			n += v[i] * v[i];
		}
		return n;
	}
	static T Norm(T* v) {
		return MathCode::sqrt(NormSquare(v));
	}
	static T Dot(T* a, T* b) {
		T d = 0;
		forN(i)
		{
			d += (a[i] * b[i]);
		}
		return d;
	}
	static void Plus(T* a, T b, T* c) {
		forN(i)
		{
			c[i] = a[i] + b;
		}
	}
	static void Minus(T* a, T b, T* c) {
		forN(i)
		{
			c[i] = a[i] - b;
		}
	}
	static void Mult(T* a, T b, T* c) {
		forN(i)
		{
			c[i] = a[i] * b;
		}
	}
	static void Divide(T* a, T b, T* c) {
		forN(i)
		{
			if (b != 0) {
				c[i] = a[i] / b;
			}
		}
	}
	void zero() {
		for (int i = 0; i < N; i++) {
			_V[i] = 0;
		}
	}
	inline T& operator[](const int j) {
		return _V[j];
	}
	inline T operator[](const int j) const {
		return _V[j];
	}
	inline T& x() {
		return _V[0];
	}
	inline T& y() {
		return _V[1];
	}
	inline T& z() {
		return _V[2];
	}
	inline T x() const {
		return _V[0];
	}
	inline T y() const {
		return _V[1];
	}
	inline T z() const {
		return _V[2];
	}
	DVec<T, N> operator+(const DVec<T, N>& a) const {
		DVec<T, N> o;
		for (int i = 0; i < N; i++) {
			o._V[i] = _V[i] + a._V[i];
		}
		return o;
	}
	DVec<T, N> operator-(const DVec<T, N>& a) const {
		DVec<T, N> o;
		for (int i = 0; i < N; i++) {
			o._V[i] = _V[i] - a._V[i];
		}
		return o;
	}
	DVec<T, N> operator*(const T a) const {
		DVec<T, N> o;
		for (int i = 0; i < N; i++) {
			o._V[i] = a * _V[i];
		}
		return o;
	}
	DVec<T, N> operator*(const DVec<T, N>& a) const {
		DVec<T, N> o;
		for (int i = 0; i < N; i++) {
			o._V[i] = _V[i] * a._V[i];
		}
		return o;
	}
	DVec<T, N> operator/(const T b) const {
		DVec<T, N> o;
		if (b == 0) {
			// treat as division by 1
			return (*this);
		}
		for (int i = 0; i < N; i++) {
			o._V[i] = _V[i] / b;
		}
		return o;
	}
	inline DVec<T, N>& proj(DVec<T, N>& v) {
		T uu = 0;
		T uv = 0;
		for (int i = 0; i < N; i++) {
			uu += (_V[i] * _V[i]);
			uv += (_V[i] * v._V[i]);
		}
		T uv_uu = uv / uu;
		if (uu > 0) {
			for (int i = 0; i < N; i++) {
				_V[i] = uv / uu * _V[i];
			}
		}
		return (*this);
	}
	DVec<T, N>& proj(DVec<T, N>& v, DVec<T, N>& r) {
		T uu = 0;
		T uv = 0;
		for (int i = 0; i < N; i++) {
			uu += (_V[i] * _V[i]);
			uv += (_V[i] * v._V[i]);
		}
		T uv_uu = uv / uu;
		if (uu > 0) {
			for (int i = 0; i < N; i++) {
				r._V[i] = uv_uu * _V[i];
			}
		}
		return (*this);
	}
	T dot(DVec<T, N>& v) {
		T innerProd = 0;
		for (int i = 0; i < N; i++) {
			innerProd += (_V[i] * v._V[i]);
		}
		return innerProd;

	}
	inline DVec<T, N>& mult(DVec<T, N>& o) {
		for (int i = 0; i < N; i++) {
			_V[i] *= o._V[i];
		}
		return (*this);
	}
	inline DVec<T, N>& mult(T a, DVec<T, N>& o) {
		for (int i = 0; i < N; i++) {
			o._V[i] = a * _V[i];
		}
		return (*this);
	}
	inline DVec<T, N>& mult(DVec<T, N>& r, DVec<T, N>& o) {
		for (int i = 0; i < N; i++) {
			o._V[i] = _V[i] * r._V[i];
		}
		return (*this);
	}
	inline DVec<T, N>& div(DVec<T, N>& o) {
		for (int i = 0; i < N; i++) {
			if (o._V[i] != 0) {
				_V[i] /= o._V[i];
			}
		}
		return (*this);
	}
	inline DVec<T, N>& div(DVec<T, N>& r, DVec<T, N>& o) {
		for (int i = 0; i < N; i++) {
			if (r._V[i] != 0) {
				o._V[i] = _V[i] / r._V[i];
			}
		}
		return (*this);
	}
	inline DVec<T, N>& plus(DVec<T, N>& o) {
		for (int i = 0; i < N; i++) {
			_V[i] += o._V[i];
		}
		return (*this);
	}
	inline DVec<T, N>& plus(DVec<T, N>& r, DVec<T, N>& o) {
		for (int i = 0; i < N; i++) {
			o._V[i] = _V[i] + r._V[i];
		}
		return (*this);
	}
	inline DVec<T, N>& minus(DVec<T, N>& o) {
		for (int i = 0; i < N; i++) {
			_V[i] -= o._V[i];
		}
		return (*this);
	}
	inline DVec<T, N>& minus(DVec<T, N>& r, DVec<T, N>& o) {
		for (int i = 0; i < N; i++) {
			o._V[i] = _V[i] - r._V[i];
		}
		return (*this);
	}

	DVec<T, N> round() const {
		DVec<T, N> o;
		for (int i = 0; i < N; i++) {
			o._V[i] = MathCode::round(_V[i]);
		}
		return o;
	}
	DVec<T, N> pow2() const {
		DVec<T, N> o;
		for (int i = 0; i < N; i++) {
			o._V[i] = _V[i] * _V[i];
		}
		return o;
	}
	T sum() const {
		T s = 0;
		for (int i = 0; i < N; i++) {
			s += _V[i];
		}
		return s;
	}
	T mag() const {
		T s = 0;
		for (int i = 0; i < N; i++) {
			s += _V[i] * _V[i];
		}
		return sqrt(s);
	}
	DVec<T, N> normal() const {
		T m = mag();
		DVec<T, N> n;
		for (int i = 0; i < N; i++) {
			n._V[i] = (_V[i] / m);
		}
		return n;
	}
	DVec<T, N> neighbor(int dir, T dist) const {
		DVec<T, N> o;
		memcpy(o._V, _V, sizeof(_V));
		o[dir] += dist;
		return o;
	}
	DVec<T, N>& set(T a) {
		_V[0] = a;
		return (*this);
	}
	DVec<T, N>& set(T a, T b) {
		_V[0] = a;
		_V[1] = b;
		return (*this);
	}
	DVec<T, N>& point(T a, T b) {
		_V[0] = a;
		_V[1] = b;
		if (N > 2) {
			_V[2] = 1;
		}
		return (*this);
	}
	DVec<T, N>& set(T a, T b, T c) {
		_V[0] = a;
		_V[1] = b;
		_V[2] = c;
		return (*this);
	}
	DVec<T, N>& point(T a, T b, T c) {
		_V[0] = a;
		_V[1] = b;
		_V[2] = c;
		if (N > 3) {
			_V[3] = 1;
		}
		return (*this);
	}
	DVec<T, N>& set(T a, T b, T c, T d) {
		_V[0] = a;
		_V[1] = b;
		_V[2] = c;
		_V[3] = d;
		return (*this);
	}
	DVec<T, N>& setN(T* a) {
		for (int i = 0; i < N; i++) {
			_V[i] = a[i];
		}
		return (*this);
	}
	DVec<T, N>& setN(T* a, int n) {
		for (int i = 0; i < n; i++) {
			_V[i] = a[i];
		}
		return (*this);
	}
	inline void copyFrom(const DVec<T, N>& o) {
		for (int i = 0; i < N; i++) {
			_V[i] = o._V[i];
		}
	}
	inline void copyTo(const DVec<T, N>& o) {
		for (int i = 0; i < N; i++) {
			o._V[i] = _V[i];
		}
	}
	inline void copyFrom(T* x) {
		for (int i = 0; i < N; i++) {
			_V[i] = x[i];
		}
	}
	inline void copyTo(T* x) {
		for (int i = 0; i < N; i++) {
			x[i] = _V[i];
		}
	}
}
;

template<typename T, int N>
class CVec: public DVec<T, N> {
private:
	T _D[N];
public:
	typedef DVec<T, N> Super;
	CVec() :
			DVec<T, N>(&_D[0]) {
		Super::zero();
	}
	virtual ~CVec() {

	}
	virtual void free() {
	}
};

typedef DVec<float, 2> Vec2;
typedef DVec<float, 3> Vec3;
typedef DVec<float, 4> Vec4;

template<typename T, int N>
inline DVec<T, N> operator*(const float a, const DVec<T, N>& v) {
	return (v * a);
}

template<typename T, int N>
ostream& operator<<(ostream& o, DVec<T, N> v) {
	o << "(" << v[0];
	for (int i = 1; i < N; i++) {
		o << "," << v[i];
	}
	o << ")";
	return o;
}

static const float _eps = 1e-6f;

template<typename T, int N>
class SquareMatrixR {
public:
	T _d[N * N];
	T* _r[N];

	SquareMatrixR() {
		for (int i = 0; i < N; i++) {
			_r[i] = &_d[i * N];
		}
		identity();
	}
	static SquareMatrixR<T, N>* New() {
		return new SquareMatrixR<T, N>();
	}
	static void Free(SquareMatrixR<T, N>* &A) {
		delete A;
		A = NULL;
	}
	inline const int size() const {
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
	inline SquareMatrixR<T, N>& setRow(int r, T* d) {
		forN(i)
		{
			_r[r][i] = d[i];
		}
		return (*this);
	}
	inline SquareMatrixR<T, N>& fillR(T* d) {
		for (int i = 0; i < size(); i++) {
			_d[i] = d[i];
		}
		return (*this);
	}
	inline SquareMatrixR<T, N>& fillC(T* d) {
		for (int j = 0; j < N; j++) {
			for (int i = 0; i < N; i++) {
				_r[j][i] = d[j * N + i];
			}
		}
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
		}
		return o;
	}
	SquareMatrixR<T, N>& inverse(SquareMatrixR<T, N>& o) {
		T det;
		inverse(o, &det);
		o.divide(det);
		return (o);
	}

	void suppressEpsilion(T eps = _eps) {
		for (int i = 0; i < N * N; i++) {
			_d[i] = (_d[i] < eps && _d[i] > -eps) ? 0.0 : _d[i];
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

	template<typename U, int M>
	SquareMatrixR<T, N>& scale(DVec<U, M>& scale) {
		if (M < N) {
			for (int i = 0; i < M; i++) {
				_r[i][i] = scale[i] * _r[i][i];
			}
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
	DVec<U, N>& mult(T* b, DVec<U, N> &c) {
		for (int i = 0; i < N; i++) {
			T e = 0;
			for (int l = 0; l < N; l++) {
				e += _r[i][l] * b[l];
			}
			c[i] = e;
		}
		return c;
	}

	DVec<T, N>& mult(DVec<T, N>& b, DVec<T, N>& c) {
		for (int i = 0; i < N; i++) {
			T e = 0;
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

	SquareMatrixR<T, N>& multC(SquareMatrixR<T, N>& b, SquareMatrixR<T, N>& c) {
		c.zero();
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				for (int k = 0; k < N; k++) {
					c._r[k][i] += _r[k][j] * b[j][i];
				}
			}
		}
		return c;
	}

	void copyFrom(SquareMatrixR<T, N>& o) {
		for (int i = 0; i < N * N; i++) {
			_d[i] = o._d[i];
		}
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

typedef SquareMatrixR<double, 2> Mat2;
typedef SquareMatrixR<double, 3> Mat3;
typedef SquareMatrixR<double, 4> Mat4;

Vec4 min(Vec4* a, int n);
Vec3 min(Vec3* a, int n);
Vec2 min(Vec2* a, int n);
Vec4 max(Vec4* a, int n);
Vec3 max(Vec3* a, int n);
Vec2 max(Vec2* a, int n);

template<typename T, int N>
DVec<T, N>& minN(DVec<T, N>* a, int n, DVec<T, N>& o) {
	o = a[0];
	for (int i = 1; i < n; i++) {
		for (int j = 0; j < N; j++) {
			o._V[j] = o._V[j] < a[i]._V[j] ? a[i]._V[j] : o._V[j];
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
	c[0] = a._d[0] * b[0] + a._d[1] * b[1] + a._d[2] * b[2];
	c[1] = a._d[3] * b[0] + a._d[4] * b[1] + a._d[5] * b[2];
	c[2] = a._d[6] * b[0] + a._d[7] * b[1] + a._d[8] * b[2];
	return c;
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

inline Mat2& mult(Mat2& a, Mat2& b, Mat2& c) {
	float e;
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

inline Mat3& mult(Mat3& a, Mat3& b, Mat3& c) {
	float e;
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

inline Mat4& mult(Mat4& a, Mat4& b, Mat4& c) {
	float e = 0;
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
