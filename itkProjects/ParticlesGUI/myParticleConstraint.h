//
//  myParticleConstraint.h
//  ParticlesGUI
//
//  Created by Joohwi Lee on 11/23/12.
//
//

#ifndef __ParticlesGUI__myParticleConstraint__
#define __ParticlesGUI__myParticleConstraint__

#include <iostream>
#include "myImageDef.h"
#include "itkSignedDanielssonDistanceMapImageFilter.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkGradientImageFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkCovariantVector.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "myParticleCore.h"

namespace pi {

    class ParticleConstraint {
    public:
        typedef itk::SignedDanielssonDistanceMapImageFilter<LabelImage, DoubleImage> SignedDistanceMapFilterType;
        typedef itk::DanielssonDistanceMapImageFilter<LabelImage, DoubleImage> DistanceMapFilterType;
        typedef DistanceMapFilterType::VectorImageType DistanceVectorImageType;
        typedef DistanceVectorImageType::PixelType DistanceVectorType;
        typedef itk::NearestNeighborInterpolateImageFunction<DoubleImage> InterpolatorType;
        typedef InterpolatorType::ContinuousIndexType ContinuousIndexType;
        typedef itk::GradientRecursiveGaussianImageFilter<LabelImage, VectorImage> GradientFilterType;
        typedef GradientFilterType::OutputPixelType GradientPixelType;
        typedef itk::VectorLinearInterpolateImageFunction<VectorImage> GradientInterpolatorType;
        typedef DistanceVectorImageType::PixelType OffsetType;

        void SetImageList(LabelVector& labelImages);

        bool IsInsideRegion(int subjId, LinearImageInterpolatorType::ContinuousIndexType& idx) const;
        bool IsInsideRegion(int subjId, LinearImageInterpolatorType::IndexType& idx) const;
        double GetDistance(int subjId, LinearImageInterpolatorType::ContinuousIndexType& idx) const;
        DistanceVectorImageType::PixelType GetInsideOffset(int subjId, DoubleImage::IndexType& idx) const;
        DistanceVectorImageType::PixelType GetOutsideOffset(int subjId, DoubleImage::IndexType& idx) const;
        GradientPixelType GetGradient(int subjId, GradientInterpolatorType::ContinuousIndexType& idx) const;

        //void ApplyConstraint(OptimizerParametersType& params) const;

        inline int GetNumberOfSubjects() const { return m_DistanceMaps.size(); }

        void Clear();

        void ApplyConstraint(ParticleSubjectArray& shapes);
        
    private:
        std::vector<DoubleImage::Pointer> m_DistanceMaps;
        std::vector<LinearImageInterpolatorType::Pointer> m_DistanceMapInterpolators;
        std::vector<DistanceVectorImageType::Pointer> m_InsideDistanceVectorMaps;
        std::vector<DistanceVectorImageType::Pointer> m_OutsideDistanceVectorMaps;
        std::vector<GradientFilterType::OutputImagePointer> m_GradientMaps;
        std::vector<GradientInterpolatorType::Pointer> m_GradientInterpolators;
    };
}
#endif /* defined(__ParticlesGUI__myParticleConstraint__) */
