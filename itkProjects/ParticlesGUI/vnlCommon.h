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
typedef vnl_matrix_fixed<double, 2, 3> VNLMat2x3;
typedef vnl_matrix_fixed<double, 2, 2> VNLMat2x2;
typedef vnl_matrix_fixed<double, 4, 4> VNLMat4x4;
typedef vnl_vector_fixed_ref<double, 2> VNLVec2Ref;
typedef vnl_vector_fixed_ref<double, 3> VNLVec3Ref;
typedef vnl_c_vector<double> VNLCVector;
typedef vnl_sse<double> VNLSSE;

typedef std::vector<VNLMatrix> VNLMatrixArray;
typedef std::vector<VNLVector> VNLVectorArray;
typedef std::vector<VNLVectorRef> VNLVectorRefArray;
typedef std::vector<double> STDDoubleArray;

class VNLMatrixRef: public VNLMatrix {
public:
    typedef VNLMatrix Base;

    // empty constructor
    VNLMatrixRef();
    VNLMatrixRef(unsigned int m, unsigned int n, double *datablck);
    VNLMatrixRef(VNLMatrixRef const & other);
    ~VNLMatrixRef();
    VNLMatrixRef& non_const() { return *this; }
    
    void reinit(unsigned int m, unsigned int n, double* datablock);
private:
    //: Resizing is disallowed
    bool resize (unsigned int, unsigned int) { return false; }
    bool make_size (unsigned int, unsigned int) { return false; }
    bool set_size (unsigned int, unsigned int) { return false; }

    //: Copy constructor from vnl_matrix<T> is disallowed
    // (because it would create a non-const alias to the matrix)
    VNLMatrixRef(VNLMatrix const &) {}
};

typedef std::vector<VNLMatrixRef> VNLMatrixRefArray;


void vnl_center(VNLMatrix& A);


//// set matrix A to identity matrix
//template <class M>
//inline void vnl_identity(M& A) {
//    int n = A.rows() > A.cols() ? A.cols() : A.rows();
//    A.fill(0);
//    for (int i = 0; i < n; i++) {
//        A[i][i] = 1;
//    }
//}

// subtract b from each row vector of A and store into X
template <class M1, class V, class M2>
void vnl_row_subtract(const M1& A, const V& b, M2& X) {
    for (int i = 0; i < A.rows(); i++) {
        for (int k = 0; k < A.cols(); k++) {
            X[i][k] = A[i][k] - b[k];
        }
    }
}

// compute row mean of A and store into b
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


// check if v has nans
bool vnl_has_nans(const VNLVector& v);


// multiply A and b and store into c
template<class M, class V1, class V2>
inline void vnl_matrix_x_vector(const M& A, const V1& b, V2& c) {
    VNLSSE::matrix_x_vector(A.begin(), b.begin(), c.begin(), A.rows(), A.cols());
}

template <class M, class V1, class V2>
inline void vnl_transform_points_2d(const M& A, const V1& b, V2& c) {
    for (int i = 0; i < b.rows(); i++) {
        c[i][0] = A[0][0] * b[i][0] + A[0][1] * b[i][1] + A[0][2];
        c[i][1] = A[1][0] * b[i][0] + A[1][1] * b[i][1] + A[1][2];
    }
}

template <class M, class V1, class V2>
inline void vnl_transform_points_3d(const M& A, const V1& b, V2& c) {
    for (int i = 0; i < b.rows(); i++) {
        c[i][0] = A[0][0] * b[i][0] + A[0][1] * b[i][1] + A[0][2] * b[i][2] + A[0][3];
        c[i][1] = A[1][0] * b[i][0] + A[1][1] * b[i][1] + A[1][2] * b[i][2] + A[1][3];
        c[i][2] = A[2][0] * b[i][0] + A[2][1] * b[i][1] + A[2][2] * b[i][2] + A[2][3];
    }
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
