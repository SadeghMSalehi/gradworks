//
//  vnlCommon.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 11/25/12.
//
//

#include "vnlCommon.h"

bool vnl_row_mean(const VNLMatrix& A, VNLVector& b)
{
    if (A.cols() != b.size()) {
        return false;
    }
    b.fill(0);
    for (int i = 0; i < A.rows(); i++) {
        for (int j = 0; j < A.cols(); j++) {
            b[j] += A[i][j];
        }
    }
    b /= A.rows();
    return true;
}

void vnl_center(VNLMatrix& A) {
    VNLVector b(A.cols());
    b.fill(0);
    for (int i = 0; i < A.rows(); i++) {
        for (int j = 0; j < A.cols(); j++) {
            b[j] += A[i][j];
        }
    }
    b /= A.rows();
    for (int i = 0; i < A.rows(); i++) {
        for (int j = 0; j < A.cols(); j++) {
            A[i][j] -= b[j];
        }
    }
    return;
}


