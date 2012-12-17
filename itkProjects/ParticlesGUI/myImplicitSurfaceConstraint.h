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
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkGradientImageFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkCovariantVector.h"
#include "itkVectorLinearInterpolateImageFunction.h"

namespace my {

class ImplicitSurfaceConstraint {
public:
    typedef itk::SignedDanielssonDistanceMapImageFilter<LabelSliceType, SliceType> SignedDistanceMapFilterType;
    typedef itk::DanielssonDistanceMapImageFilter<LabelSliceType, SliceType> DistanceMapFilterType;
    typedef DistanceMapFilterType::VectorImageType DistanceVectorImageType;
    typedef DistanceVectorImageType::PixelType DistanceVectorType;
    typedef itk::NearestNeighborInterpolateImageFunction<SliceType> InterpolatorType;
    typedef InterpolatorType::ContinuousIndexType ContinuousIndexType;
    typedef itk::CovariantVector<double,2> GradientVectorType;
    typedef itk::Image<GradientVectorType,2> GradientImageType;
    typedef itk::GradientRecursiveGaussianImageFilter<LabelSliceType, GradientImageType> GradientFilterType;
    typedef GradientFilterType::OutputPixelType GradientPixelType;
    typedef itk::VectorLinearInterpolateImageFunction<GradientImageType> GradientInterpolatorType;
    typedef DistanceVectorImageType::PixelType OffsetType;

    void SetImageList(ImageContainer::List* imageList);

    bool IsInsideRegion(int subjId, SliceInterpolatorType::ContinuousIndexType& idx) const;
    bool IsInsideRegion(int subjId, SliceInterpolatorType::IndexType& idx) const;
    double GetDistance(int subjId, SliceInterpolatorType::ContinuousIndexType& idx) const;
    DistanceVectorImageType::PixelType GetInsideOffset(int subjId, SliceType::IndexType& idx) const;
    DistanceVectorImageType::PixelType GetOutsideOffset(int subjId, SliceType::IndexType& idx) const;
    GradientPixelType GetGradient(int subjId, GradientInterpolatorType::ContinuousIndexType& idx) const;

    void ApplyConstraint(OptimizerParametersType& params) const;

    inline int GetNumberOfSubjects() const { return m_DistanceMaps.size(); }

    void Clear();
private:
    std::vector<SliceType::Pointer> m_DistanceMaps;
    std::vector<InterpolatorType::Pointer> m_DistanceMapInterpolators;
    std::vector<DistanceVectorImageType::Pointer> m_InsideDistanceVectorMaps;
    std::vector<DistanceVectorImageType::Pointer> m_OutsideDistanceVectorMaps;
    std::vector<GradientFilterType::OutputImagePointer> m_GradientMaps;
    std::vector<GradientInterpolatorType::Pointer> m_GradientInterpolators;
};
}
#endif /* defined(__ParticlesGUI__myImplicitSurfaceConstraint__) */
