//
//  ExpLogTransform.h
//  pnsc++
//
//  Created by Joohwi Lee on 10/12/12.
//
//

#ifndef __pnsc____ExpLogTransform__
#define __pnsc____ExpLogTransform__

#include <iostream>
#include "PNSBase.h"
#include "iostream"

class ExpLogTransform {
public:
    typedef PNSBase::VectorType VectorType;
    typedef PNSBase::MatrixType MatrixType;

    // Tangent space to sphere
    void TransformExp(MatrixType& t, MatrixType& s);

    // Sphere to tangent space
    void TransformLog(MatrixType& s, MatrixType& t);

    // Set tangent point
    void SetTangentPoint(VectorType& p);
    VectorType& GetTangentPoint();
    MatrixType& GetRotationMatrix();


private:    
    VectorType m_TangentPoint;
    MatrixType m_Rotation;

};
#endif /* defined(__pnsc____ExpLogTransform__) */
