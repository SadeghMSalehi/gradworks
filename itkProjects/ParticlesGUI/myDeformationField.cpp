//
//  myDeformationField.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 12/6/12.
//
//

#include "myDeformationField.h"
#include "itkImageRegionConstIteratorWithIndex.h"

myDeformationField::myDeformationField() {

}

myDeformationField::~myDeformationField() {

}

void myDeformationField::SetLandmarks(int n, double* src, double* dst) {
    m_nPoints = n;
    m_Src.reinit(n, SDim, src);
    m_Dst.reinit(n, SDim, dst);
    m_Displacement = m_Dst - m_Src;
}

SliceType::Pointer myDeformationField::ResampleLinear(SliceType::Pointer srcImg, LabelSliceType::Pointer roiImg) {
    LabelSliceIteratorType labelIter(roiImg, roiImg->GetBufferedRegion());
    for (labelIter.GoToBegin(); labelIter.IsAtEnd(); ++labelIter) {
        // if the voxel is ROI
        if (labelIter.Get() > 0) {
            VNLVec2 disp;
            disp.fill(0);
            for (int i = 0; i < m_nPoints; i++) {
                
            }
        }
    }
    return SliceType::Pointer(NULL);
}
