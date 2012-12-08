//
//  myBSplineRegistration.h
//  ParticlesGUI
//
//  Created by Joohwi Lee on 12/6/12.
//
//

#ifndef __ParticlesGUI__myBSplineRegistration__
#define __ParticlesGUI__myBSplineRegistration__

#include <iostream>
#include "vnlCommon.h"
#include "itkPointSet.h"
#include "itkBSplineScatteredDataPointSetToImageFilter.h"
#include "myImageContainer.h"

namespace my {
    typedef itk::BSplineScatteredDataPointSetToImageFilter
    <DisplacementFieldPointSetType, DisplacementFieldType> BSplineFilterType;
    typedef BSplineFilterType::WeightsContainerType WeightsContainerType;

    class BSplineRegistration {
    public:
        BSplineRegistration();
        ~BSplineRegistration();

        void SetReferenceImage(SliceType::Pointer refImage);
        void SetLandmarks(int n, double* src, double* dst);
        void Update();
        SliceType::Pointer WarpImage(SliceType::Pointer srcImage);
        SliceType::Pointer GetDisplacementMagnitude();
        DisplacementFieldType::Pointer GetDisplacementField();
        DisplacementFieldType::Pointer GetControlPoints();

    private:
        VNLVector m_Params;
        SliceType::Pointer m_RefImage;
        DisplacementFieldPointSetType::Pointer m_FieldPoints;
        DisplacementFieldType::Pointer m_DisplacementField;
        DisplacementFieldType::Pointer m_PhiLattice;
    };
}

#endif /* defined(__ParticlesGUI__myBSplineRegistration__) */