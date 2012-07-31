/*
 * GramSchmidt.h
 *
 *  Created on: Jul 25, 2012
 *      Author: joohwile
 */

#ifndef GRAMSCHMIDT_H_
#define GRAMSCHMIDT_H_

#include "MatrixCode.h"

namespace MathCode {
template<typename T, int N>
void row_copy(SquareMatrixR<T, N>& src, int j, SquareMatrixR<T, N>& dst) {
    //memcpy(dst, src, sizeof(T) * j);
    for (int i = 0; i < N; i++) {
        dst._r[j][i] = src._r[j][i];
    }
}

template<typename T, int N>
void gramschmidt(SquareMatrixR<T, N>& src, SquareMatrixR<T, N>& dst) {
    DVec<T, N> u0(dst._r[0]);
    u0.copyFrom(src._r[0]);
    DVec<T, N> proj;
    for (int i = 1; i < N; i++) {
        DVec<T, N> ui(dst._r[i]), vi(src._r[i]);
        ui.copyFrom(src._r[i]);
        for (int j = 0; j < i; j++) {
            DVec<T, N> uj(dst._r[j]);
            uj.proj(vi, proj);
            ui.minus(proj);
        }
    }
}

template<typename T, int N>
void gramschmidtQR(SquareMatrixR<T, N>& W, SquareMatrixR<T, N>& Q, SquareMatrixR<T, N>& R) {
    R.zero();
    R._r[0][0] = DVec<T, N>::Norm(W._r[0]);
    DVec<T, N>::Divide(W._r[0], R._r[0][0], Q._r[0]);
    for (int i = 1; i < N; i++) {
        T riri = DVec<T, N>::NormSquare(W._r[i]);
        Q.setRow(i, W._r[i]);
        for (int j = 0; j < i; j++) {
            R._r[i][j] = DVec<T, N>::Dot(W._r[i], Q._r[j]);
            riri -= (R._r[i][j] * R._r[i][j]);
            forN(k) {
                Q._r[i][k] -= R._r[i][j] * Q._r[j][k];
            }
        }
        R._r[i][i] = sqrt(riri);
        DVec<T, N>::Divide(Q._r[i], R._r[i][i], Q._r[i]);
    }
}

template<typename T, int N>
void gramschmidtQRx(SquareMatrixR<T, N>& W, SquareMatrixR<T, N>& Q, SquareMatrixR<T, N>& R) {
    SquareMatrixR<T, N> Wx;
    Wx.fillR(W._d);
    R.zero();
    DVec<T, N>::Divide(W._r[0], R._r[0][0], Q._r[0]);
    for (int i = 0; i < N; i++) {
        R._r[i][i] = DVec<T, N>::Norm(Wx._r[i]);
        DVec<T, N>::Divide(Wx._r[i], R._r[i][i], Q._r[i]);
        for (int j = i + 1; j < N; j++) {
            R._r[j][i] = DVec<T, N>::Dot(Wx._r[j], Q._r[i]);
            forN(k) {
                Wx._r[j][k] -= (R._r[j][i] * Q._r[i][k]);
            }
        }
    }
}
}

#endif /* GRAMSCHMIDT_H_ */
