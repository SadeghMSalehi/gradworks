//
//  vnlCommon.h
//  ParticlesGUI
//
//  Created by Joohwi Lee on 11/25/12.
//
//

#ifndef __ParticlesGUI__vnlCommon__
#define __ParticlesGUI__vnlCommon__

#include <iostream>
#include "vnl/vnl_vector.h"
#include "vnl/vnl_vector_fixed.h"
#include "vnl/vnl_vector_ref.h"
#include "vnl/vnl_vector_fixed_ref.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_diag_matrix.h"
#include "vnl/vnl_matrix_ref.h"
#include "vnl/vnl_c_vector.h"
#include "vnl/vnl_sse.h"

#include "vector"

typedef vnl_vector<double> VNLVector;
typedef vnl_matrix<double> VNLMatrix;
typedef vnl_diag_matrix<double> VNLDiagMatrix;
typedef vnl_vector_ref<double> VNLVectorRef;
typedef vnl_vector_fixed<double, 3> VNLVec3;
typedef vnl_vector_fixed<double, 2> VNLVec2;
typedef vnl_matrix_fixed<double, 3, 3> VNLMat3x3;
typedef vnl_matrix_fixed<double, 4, 4> VNLMat4x4;
typedef vnl_vector_fixed_ref<double, 2> VNLVec2Ref;
typedef vnl_vector_fixed_ref<double, 3> VNLVec3Ref;
typedef vnl_matrix_ref<double> VNLMatrixRef;
typedef vnl_c_vector<double> VNLCVector;
typedef vnl_sse<double> VNLSSE;

typedef std::vector<VNLMatrix> VNLMatrixArray;
typedef std::vector<VNLVector> VNLVectorArray;
typedef std::vector<VNLVectorRef> VNLVectorRefArray;
typedef std::vector<VNLMatrixRef> VNLMatrixRefArray;
typedef std::vector<double> STDDoubleArray;

void vnl_center(VNLMatrix& A);

template <class M>
inline void vnl_identity(M& A) {
    int n = A.rows() > A.cols() ? A.cols() : A.rows();
    A.fill(0);
    for (int i = 0; i < n; i++) {
        A[i][i] = 1;
    }
}

template <class M1, class V, class M2>
void vnl_row_subtract(const M1& A, const V& b, M2& X) {
    for (int i = 0; i < A.rows(); i++) {
        for (int k = 0; k < A.cols(); k++) {
            X[i][k] = A[i][k] - b[k];
        }
    }
}

template <class M, class V>
bool vnl_row_mean(const M& A, V& b) {
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

bool vnl_has_nan(const VNLVector& v);


template<class M, class V1, class V2>
inline void vnl_matrix_x_vector(const M& A, const V1& b, V2& c) {
    VNLSSE::matrix_x_vector(A.begin(), b.begin(), c.begin(), A.rows(), A.cols());
}


struct VNLAlgebra {
    int R;
    int C;
    VNLAlgebra(int r, int c) { R = r; C = c; }
    inline void multiply_matrix_vector(const double* A, const double* b, double* c) {
        VNLSSE::matrix_x_vector(A, b, c, R, C);
    }
};

#endif /* defined(__ParticlesGUI__vnlCommon__) */
