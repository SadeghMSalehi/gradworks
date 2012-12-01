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

#include "vector"

typedef vnl_vector<double> VNLVector;
typedef vnl_matrix<double> VNLMatrix;
typedef vnl_diag_matrix<double> VNLDiagMatrix;
typedef vnl_vector_ref<double> VNLVectorRef;
typedef vnl_vector_fixed<double, 3> VNLVec3;
typedef vnl_vector_fixed<double, 2> VNLVec2;
typedef vnl_vector_fixed_ref<double, 2> VNLVec2Ref;
typedef vnl_vector_fixed_ref<double, 3> VNLVec3Ref;
typedef vnl_matrix_ref<double> VNLMatrixRef;
typedef vnl_c_vector<double> VNLCVector;

typedef std::vector<VNLMatrix> VNLMatrixArray;
typedef std::vector<VNLVector> VNLVectorArray;
typedef std::vector<VNLVectorRef> VNLVectorRefArray;
typedef std::vector<VNLMatrixRef> VNLMatrixRefArray;

void vnl_center(VNLMatrix& A);
bool vnl_row_mean(const VNLMatrix& A, VNLVector& b);
bool vnl_has_nan(const VNLVector& v);



#endif /* defined(__ParticlesGUI__vnlCommon__) */
