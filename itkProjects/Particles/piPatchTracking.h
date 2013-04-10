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
#include <QPolygonF>

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
        void setPatchRegion(RealImage::RegionType region);
        void translatePatchRegion(RealImage::RegionType region);
        void resamplePatch();

        void beginTracking();
        void stepTracking();

        RealImage::Pointer getPatch(int i);
        QPolygonF& getPatchPolygon(int i);

    private:
        void extractPatch(int i, bool newPatch = true);
        void transformPatchRegion();
        void setupOptimizer(itk::RegularStepGradientDescentOptimizer* opti);
        void setupOptimizer(itk::FRPROptimizer* opti);

    private:
        RealImage::Pointer _images[2];
        RealImage::Pointer _patches[2];
        RealImage::RegionType _patchRegion[2];
        QPolygonF _patchPolygon[2];

        typedef itk::TranslationTransform<double,2> TransformType;
//        typedef itk::Rigid2DTransform<double> TransformType;
        TransformType::Pointer _transform;

        typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
        OptimizerType::Pointer _optimizer;

        typedef itk::NormalizedCorrelationImageToImageMetric<RealImage, RealImage> CostFunctionType;
        CostFunctionType::Pointer _costFunc;

        typedef OptimizerType::ParametersType ParametersType;
    };
}
#endif /* defined(__ParticleGuidedRegistration__piPatchTracking__) */