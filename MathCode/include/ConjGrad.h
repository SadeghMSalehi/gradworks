/*
 * ConjGrad.h
 *
 *  Created on: Jul 25, 2012
 *      Author: joohwile
 */

#ifndef CONJGRAD_H_
#define CONJGRAD_H_

#include "MatrixCode.h"

namespace MathCode {

template<typename T, int N>
class ConjGrad {
public:
    typedef DVec<T, N> Vec;
    typedef SquareMatrixR<T, N> Mat;

private:
    Mat* _A;
    Vec* _b, *_x0, *_xs;

public:
    ConjGrad(Mat* A, Vec* b, Vec* x0, Vec* xs) {
        _A = A;
        _b = b;
        _x0 = x0;
        _xs = xs;
    }

    T computeAk(Vec& rk, Vec& pk, Vec& Apk) {
        T pkApk = Apk.dot(pk);
        T rkrk = rk.dot(rk);
        return rkrk / pkApk;
    }

    T computeBk(Vec& rkk, Vec& rk) {
        double rkkrkk = rkk.dot(rkk);
        double rkrk = rk.dot(rk);
        return rkkrkk / rkrk;
    }

    bool checkTermination(Vec& rkk, T eps = 1e-5) {
        return rkk.mag() < eps;
    }

    void compute() {
        Vec xk(_x0), xkk, rk, rkk, pk, pkk, apk, bpk, aApk, Apk;

        Mat & A = (*_A);
        Vec& b = (*_b);

        A.mult(xk, rk);
        b.minus(rk, rk);
        pk.copyFrom(rk._V);

        for (int k = 0; k < N; k++) {
            _DBG_(cout << "xk = " << xk << "; rk = " << rk << "; pk = " << pk << ";" << endl);
            A.mult(pk, Apk);
            T ak = computeAk(rk, pk, Apk);
            pk.mult(ak, apk);
            Apk.mult(ak, aApk);
            xk.plus(apk, xkk);
            rk.minus(aApk, rkk);
            if (checkTermination(rkk)) {
                break;
            }
            T bk = computeBk(rkk, rk);
            pk.mult(bk, bpk);
            rkk.plus(bpk, pkk);
            pk.copyFrom(pkk._V);
            rk.copyFrom(rkk._V);
            xk.copyFrom(xkk._V);
        }

        _xs->copyFrom(xkk._V);
    }
};
}
#endif /* CONJGRAD_H_ */
