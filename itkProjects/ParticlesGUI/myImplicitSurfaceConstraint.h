//
//  myImplicitSurfaceConstraint.h
//  ParticlesGUI
//
//  Created by Joohwi Lee on 11/23/12.
//
//

#ifndef __ParticlesGUI__myImplicitSurfaceConstraint__
#define __ParticlesGUI__myImplicitSurfaceConstraint__

#include <iostream>
#include "itkOptimizerCommon.h"
#include "myImageContainer.h"

#include "itkSignedDanielssonDistanceMapImageFilter.h"

class myImplicitSurfaceConstraint {
public:
    typedef itk::SignedDanielssonDistanceMapImageFilter<SliceType, SliceType> DistanceMapFilterType;
    typedef DistanceMapFilterType::VectorImageType DistanceVectorImageType;

    void AddSurface(SliceType::Pointer labelMap);
    double GetDistance(int subjId, SliceInterpolatorType::ContinuousIndexType& idx);
    DistanceVectorImageType::PixelType GetDistanceVector(int subjId, SliceType::IndexType& idx);
    void ApplyConstraint(OptimizerParametersType& params);

    inline int GetNumberOfSubjects() { return m_DistanceMaps.size(); }
private:
    std::vector<SliceType::Pointer> m_DistanceMaps;
    std::vector<SliceInterpolatorType::Pointer> m_DistanceMapInterpolators;
    std::vector<DistanceVectorImageType::Pointer> m_DistanceVectorMaps;
};
#endif /* defined(__ParticlesGUI__myImplicitSurfaceConstraint__) */
