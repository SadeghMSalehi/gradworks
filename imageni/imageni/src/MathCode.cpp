/*
 * CMath.cpp
 *
 *  Created on: Jul 21, 2012
 *      Author: joohwile
 */

#include "MathCode.h"
#include "math.h"

namespace cmath {

float abs(float n) {
    return n > 0 ? n : -n;
}
int round(float n) {
    return (int) (n > 0 ? (n + .5) : (n - .5));
}
float sqrt(float s) {
    return ::sqrt(s);
}
float cos(float x) {
    return ::cos(x);
}
float sin(float x) {
    return ::sin(x);
}
Vec2 max(Vec2* a, int n) {
    Vec2 o = a[0];
    for (int i = 0; i < n; i++) {
        o._x[0] = o._x[0] < a[i]._x[0] ? a[i]._x[0] : o._x[0];
        o._x[1] = o._x[1] < a[i]._x[1] ? a[i]._x[1] : o._x[1];
    }
    return o;
}

Vec3 max(Vec3* a, int n) {
    Vec3 o = a[0];
    for (int i = 0; i < n; i++) {
        o._x[0] = o._x[0] < a[i]._x[0] ? a[i]._x[0] : o._x[0];
        o._x[1] = o._x[1] < a[i]._x[1] ? a[i]._x[1] : o._x[1];
        o._x[2] = o._x[2] < a[i]._x[2] ? a[i]._x[2] : o._x[2];
    }
    return o;
}

Vec4 max(Vec4* a, int n) {
    Vec4 o = a[0];
    for (int i = 0; i < n; i++) {
        o._x[0] = o._x[0] < a[i]._x[0] ? a[i]._x[0] : o._x[0];
        o._x[1] = o._x[1] < a[i]._x[1] ? a[i]._x[1] : o._x[1];
        o._x[2] = o._x[2] < a[i]._x[2] ? a[i]._x[2] : o._x[2];
        o._x[3] = o._x[3] < a[i]._x[3] ? a[i]._x[3] : o._x[3];
    }
    return o;
}

Vec2 min(Vec2* a, int n) {
    Vec2 o = a[0];
    for (int i = 0; i < n; i++) {
        o._x[0] = o._x[0] > a[i]._x[0] ? a[i]._x[0] : o._x[0];
        o._x[1] = o._x[1] > a[i]._x[1] ? a[i]._x[1] : o._x[1];
    }
    return o;
}

Vec3 min(Vec3* a, int n) {
    Vec3 o = a[0];
    for (int i = 0; i < n; i++) {
        o._x[0] = o._x[0] > a[i]._x[0] ? a[i]._x[0] : o._x[0];
        o._x[1] = o._x[1] > a[i]._x[1] ? a[i]._x[1] : o._x[1];
        o._x[2] = o._x[2] > a[i]._x[2] ? a[i]._x[2] : o._x[2];
    }
    return o;
}

Vec4 min(Vec4* a, int n) {
    Vec4 o = a[0];
    for (int i = 0; i < n; i++) {
        o._x[0] = o._x[0] > a[i]._x[0] ? a[i]._x[0] : o._x[0];
        o._x[1] = o._x[1] > a[i]._x[1] ? a[i]._x[1] : o._x[1];
        o._x[2] = o._x[2] > a[i]._x[2] ? a[i]._x[2] : o._x[2];
        o._x[3] = o._x[3] > a[i]._x[3] ? a[i]._x[3] : o._x[3];
    }
    return o;
}

template <typename T>
inline T dot(T* a, T* b, int n) {
	T c = 0;
	for (int i = 0; i < n; i++) {
		c += (a[i]*b[i]);
	}
	return c;
}

template <typename T>
inline T dot(MatrixR<T> &a, int row, MatrixR<T> &b, int col) {
	T c = 0;
	for (int i = 0; i < b.nRows; i++) {
		c += (a[row][i]*b[i][col]);
	}
	return c;
}

template <typename T>
void mult(MatrixR<T>& a, MatrixR<T>& b, MatrixR<T>& o) {
	for (int i = 0; i < a.nRows; i++) {
		for (int j = 0; j < b.nCols; j++) {
			T c = 0;
			for (int k = 0; k < b.nRows; k++) {
				c += (a[i][k]*b[k][j]);
			}
			o[i][j] = c;
		}
	}
}

}
