//
//  pnsMath.h
//  pnsc++
//
//  Created by Joohwi Lee on 10/11/12.
//
//

#ifndef __pnsc____pnsMath__
#define __pnsc____pnsMath__

#include <iostream>
#include "PNSBase.h"

class PNSMath {
public:
    typedef PNSBase::VectorType VectorType;
    typedef PNSBase::MatrixType MatrixType;

    MatrixType m_Data;
    VectorType m_Normal;
    VectorType m_CenterAtTangent;
    double m_Phi;

    PNSMath() {

    }
    void startOptimization(VectorType n0, double phi);
};

#endif /* defined(__pnsc____pnsMath__) */
