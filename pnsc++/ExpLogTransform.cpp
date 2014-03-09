//
//  ExpLogTransform.cpp
//  pnsc++
//
//  Created by Joohwi Lee on 10/12/12.
//
//

#include "ExpLogTransform.h"
#include "math.h"

// Tangent space to sphere
void ExpLogTransform::TransformExp(MatrixType& t, MatrixType& s) {
    int nCols = t.n_cols;
    int nDim  = t.n_rows;
    MatrixType x;
    x.zeros(nDim, nCols);
    for (int j = 0; j < nCols; j++) {
        double zMag = 0;
        for (int i = 0; i < nDim - 1; i++) {
            zMag += (t.at(i,j) * t.at(i,j));
        }
        zMag = sqrt(zMag);
        for (int i = 0; i < nDim - 1; i++) {
            x.at(i,j) = sin(zMag) / zMag * t.at(i,j);
        }
        x.at(nDim-1,j) = cos(zMag);
    }
    s = m_Rotation.i() * x;
}

// Sphere to tangent space
void ExpLogTransform::TransformLog(MatrixType& s, MatrixType& t) {
    MatrixType x = m_Rotation * s;

//     x.print("Rotated Data: ");
    int nCols = s.n_cols;
    int nDim  = s.n_rows;
    t.zeros(nDim, nCols);

    for (int j = 0; j < nCols; j++) {
        double theta = acos(x.at(nDim-1,j));
        for (int i = 0; i < nDim - 1; i++) {
            t.at(i,j) = theta / sin(theta) * x.at(i,j);
        }
        t.at(nDim-1,j) = 1; //x[nDim-1][j];
    }
}

// Set tangent point
void ExpLogTransform::SetTangentPoint(VectorType& p) {
    m_TangentPoint = p;
    m_TangentPoint /= arma::norm(p, 2);
    VectorType northPole(p.n_rows);
    northPole(p.n_rows - 1) = 1;
    PNSBase::ComputeRotationMatrix(m_TangentPoint, northPole, 0, m_Rotation);
}

ExpLogTransform::VectorType& ExpLogTransform::GetTangentPoint() {
    return m_TangentPoint;
}

ExpLogTransform::MatrixType& ExpLogTransform::GetRotationMatrix() {
    return m_Rotation;
}