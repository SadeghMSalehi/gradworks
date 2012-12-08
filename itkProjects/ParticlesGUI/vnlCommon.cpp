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

    // debug: check row mean b is correct
    //    std::cout << "Mean: " << b << std::endl;
    return;
}

bool vnl_has_nans(const VNLVector& params) {
    for (int i = 0; i < params.size(); i++) {
        if (params[i] != params[i]) {
            return true;
        }
    }
    return false;
}


VNLMatrixRef::VNLMatrixRef() : VNLMatrix() {

}

VNLMatrixRef::VNLMatrixRef(unsigned int m, unsigned int n, double *datablck) {
    Base::data = VNLCVector::allocate_Tptr(m);
    for (unsigned int i = 0; i < m; ++i)
        Base::data[i] = datablck + i * n;
    Base::num_rows = m;
    Base::num_cols = n;
#if VCL_HAS_SLICED_DESTRUCTOR_BUG
    this->vnl_matrix_own_data = 0;
#endif
}

VNLMatrixRef::VNLMatrixRef(VNLMatrixRef const & other) : VNLMatrix() {
    Base::data = VNLCVector::allocate_Tptr(other.rows());
    for (unsigned int i = 0; i < other.rows(); ++i)
        Base::data[i] = const_cast<double*>(other.data_block()) + i * other.cols();
    Base::num_rows = other.rows();
    Base::num_cols = other.cols();
#if VCL_HAS_SLICED_DESTRUCTOR_BUG
    this->vnl_matrix_own_data = 0;
#endif
}

VNLMatrixRef::~VNLMatrixRef() {
    Base::data[0] = 0; // Prevent base dtor from releasing our memory
}

void VNLMatrixRef::reinit(unsigned int m, unsigned int n, double* datablck) {
    if (Base::num_rows != m) {
        if (Base::data != NULL) {
            VNLCVector::deallocate(Base::data, Base::num_rows);
            Base::data = VNLCVector::allocate_Tptr(m);
        }
    }
    for (unsigned int i = 0; i < m; ++i) {
        Base::data[i] = datablck + i * n;
    }
    Base::num_rows = m;
    Base::num_cols = n;
#if VCL_HAS_SLICED_DESTRUCTOR_BUG
    this->vnl_matrix_own_data = 0;
#endif
}