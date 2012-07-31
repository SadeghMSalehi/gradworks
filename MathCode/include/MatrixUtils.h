/*
 * MatrixUtils.h
 *
 *  Created on: Jul 29, 2012
 *      Author: joohwile
 */

#ifndef MATRIXUTILS_H_
#define MATRIXUTILS_H_

#include "MatrixCode.h"


namespace MathCode {


template<typename T, int N>
class MatrixFunctions {
public:
    typedef SquareMatrixR<T, N> Matrix;
    typedef DVec<T, N> Vector;

    static void diagonal(Matrix* A, Matrix* D) {
        forN(i) {
            D->_r[i][i] = A->_r[i][i];
        }
    }
    static void off_diagonal(Matrix* A, Matrix* B) {
        forN(j) {
            forN(i) {
                if (i != j) {
                    B->_r[i][j] = A->_r[i][j];
                } else {
                    B->_r[i][j] = 0;
                }
            }
        }
    }

    static void random(Matrix* A) {
        for (int i = 0; i < N * N; i++) {
            A->_d[i] = ::rand() / 1e9;
        }
    }
    static void image2standard(Matrix* A, int w, int h) {
        if (N == 3) {
            //_d[0] = -_d[0];
            //_d[2] = _d[2] + w;
            A->_d[4] = -A->_d[4];
            A->_d[5] = A->_d[5] + h;
        }
    }
    static void standard2image(Matrix* A, int w, int h) {
        if (N == 3) {
            A->_d[4] = -A->_d[4];
            A->_d[5] = A->_d[5] - h;
        }
    }
    void rotateAxisAround(Matrix* A, T deg, T cx, T cy) {
        T eps = 1;
        if (deg < eps && deg > eps) {
            return;
        }
        T rad = deg / 180 * PI;
        T cr = cos(rad), sr = sin(rad);
        if (N == 3) {
            T a = A->_d[0];
            T b = A->_d[1];
            T c = A->_d[2] - cx;
            T d = A->_d[3];
            T e = A->_d[4];
            T f = A->_d[5] - cy;
            A->_d[0] = a * cr + d * sr;
            A->_d[1] = b * cr + e * sr;
            A->_d[2] = c * cr + f * sr + cx;
            A->_d[3] = a * -sr + d * cr;
            A->_d[4] = b * -sr + e * cr;
            A->_d[5] = c * -sr + f * cr + cy;
        } else if (N == 2) {
            T a = A->_d[0];
            T b = A->_d[1];
            T c = A->_d[2];
            T d = A->_d[3];
            A->_d[0] = cr * a + sr * c;
            A->_d[1] = cr * b + sr * d;
            A->_d[2] = -sr * a + cr * c;
            A->_d[3] = -sr * b + cr * d;
        }
    }
    void translateOrigin(Matrix* A, Vec2& v) {
        A->_r[0][N - 1] = A->_r[0][N - 1] - v[0];
        A->_r[1][N - 1] = A->_r[1][N - 1] - v[1];
    }
    void translateOrigin(Matrix* A, double tx, double ty) {
        A->_r[0][N - 1] = A->_r[0][N - 1] - tx;
        A->_r[1][N - 1] = A->_r[1][N - 1] - ty;
    }
    void translateOrigin(Matrix* A, Vec3& v) {
        A->_r[0][N - 1] = A->_r[0][N - 1] - v[0];
        A->_r[1][N - 1] = A->_r[1][N - 1] - v[1];
        A->_r[2][N - 1] = A->_r[2][N - 1] - v[2];
    }
    void translateOriginN(Matrix* A, DVec<T, N>& v) {
        for (int i = 0; i < N - 1; i++) {
            A->_r[i][N - 1] = A->_r[i][N - 1] - v[i];
        }
    }
};

}
#endif /* MATRIXUTILS_H_ */
