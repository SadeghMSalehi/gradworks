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
#include "vnl/vnl_matrix.h"
#include "vector"

typedef vnl_vector<double> VNLVector;
typedef vnl_matrix<double> VNLMatrix;
typedef std::vector<VNLMatrix> VNLMatrixArray;
typedef std::vector<VNLVector> VNLVectorRefArray;

void vnl_center(VNLMatrix& A);
bool vnl_row_mean(const VNLMatrix& A, VNLVector& b);

#endif /* defined(__ParticlesGUI__vnlCommon__) */
