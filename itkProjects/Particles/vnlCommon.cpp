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

void vnl_row_center(VNLMatrix& A) {
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


// compute covariance matrix of A
// A is row-major matrix
void vnl_covariance(const VNLMatrix& A, VNLMatrix& cov) {
    const int nCols = A.cols();
    cov.set_size(A.rows(), A.rows());
    for (int i = 0; i < cov.rows(); i++) {
        for (int j = i; j < cov.cols(); j++) {
            cov[j][i] = cov[i][j] = VNLCVector::dot_product(A[i], A[j], nCols);
        }
    }
}

void vnl_add_diagonal(VNLMatrix& A, double alpha) {
    const int n = (A.rows() < A.cols()) ? A.rows() : A.cols();
    for (int i = 0; i < n; i++) {
        A[i][i] += n;
    }
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