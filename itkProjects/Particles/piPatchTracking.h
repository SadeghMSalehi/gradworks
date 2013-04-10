//
//  piPatchTracking.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/9/13.
//
//

#ifndef __ParticleGuidedRegistration__piPatchTracking__
#define __ParticleGuidedRegistration__piPatchTracking__

#include <iostream>
#include <itkCenteredRigid2DTransform.h>
#include <itkTranslationTransform.h>
#include <itkNormalizedCorrelationImageToImageMetric.h>
#include <itkMeanSquaresImageToImageMetric.h>
#include <itkRegularStepGradientDescentOptimizer.h>
#include <itkLBFGSBOptimizer.h>
#include <itkFRPROptimizer.h>

#include "piImageDef.h"

namespace pi {
    class PatchTracking {
    public:
        PatchTracking();
        ~PatchTracking() {}

        inline void operator()() {
            stepTracking();
        }

        void setImage(int i, RealImage::Pointer image);
        void setInitialRegion(int i, RealImage::RegionType region);

        void beginTracking();
        void stepTracking();

        RealImage::Pointer getPatch(int i);
        RealImage::IndexType& getFinalIndex(int i);

    private:
        void extractPatch(int i);

    private:
        RealImage::Pointer _images[2];
        RealImage::Pointer _patches[2];
        RealImage::RegionType _initialRegion[2];
        RealImage::PointType _centerOfRegion[2];
        RealImage::PointType _finalPoints[2];
        RealImage::IndexType _finalIndex[2];

        typedef itk::TranslationTransform<double,2> TransformType;
        TransformType::Pointer _transform;

        typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
        OptimizerType::Pointer _optimizer;

        typedef itk::MeanSquaresImageToImageMetric<RealImage, RealImage> CostFunctionType;
        CostFunctionType::Pointer _costFunc;

        typedef OptimizerType::ParametersType ParametersType;
    };
}
#endif /* defined(__ParticleGuidedRegistration__piPatchTracking__) */