//
//  myDeformationField.h
//  ParticlesGUI
//
//  Created by Joohwi Lee on 12/6/12.
//
//

#ifndef __ParticlesGUI__myDeformationField__
#define __ParticlesGUI__myDeformationField__

#include <iostream>
#include "myImageContainer.h"
#include "vnlCommon.h"

class myDeformationField {
public:
    myDeformationField();
    ~myDeformationField();

    void SetLandmarks(int n, double* src, double* dst);
    SliceType::Pointer ResampleLinear(SliceType::Pointer srcImg, LabelSliceType::Pointer label);

private:
    int m_nPoints;
    VNLMatrixRef m_Src;
    VNLMatrixRef m_Dst;
    VNLMatrix m_Displacement;
};

#endif /* defined(__ParticlesGUI__myDeformationField__) */
