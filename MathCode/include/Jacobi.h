/*
 * Jacobi.h
 *
 *  Created on: Jul 29, 2012
 *      Author: joohwile
 */

#ifndef JACOBI_H_
#define JACOBI_H_

#include "MatrixCode.h"
#include "MatrixUtils.h"

namespace MathCode {

template<typename T, int N>
void jacobi_solve(SquareMatrixR<T, N>* A, DVec<T, N>* b, DVec<T, N>* x0, DVec<T, N>* xs) {
    T eps = 1e-5;
    T maxError = MAX_T;
    _DBG_(cout << "A = " << *A << endl);
    _DBG_(cout << "x0 = " << *x0 << endl);
    for (int iter = 0; maxError > eps && iter < 10; iter++) {
        forN(j) {
            T d = 0;
            if (A->_r[j][j] != 0) {
                d = 1.f / A->_r[j][j];
            }
            xs->_V[j] = 0;
            forN(i) {
                if (i != j) {
                    xs->_V[j] += (-d * A->_r[j][i] * x0->_V[i]);
                }
            }
            xs->_V[j] += d * b->_V[j];
            maxError = __mathcode_max(maxError, abs(xs->_V[j]-x0->_V[j]));
        }
        _DBG_(cout << "Iteration: " << iter << "; xs = " << *xs << endl);
        x0->copyFrom(xs->_V);
    }
}

}

#endif /* JACOBI_H_ */
