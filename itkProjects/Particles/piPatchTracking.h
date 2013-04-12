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
#include <itkBSplineTransform.h>
#include <itkMeanSquaresImageToImageMetricv4.h>
#include <itkMattesMutualInformationImageToImageMetricv4.h>
#include "itkEntropyImageToImageMetricv4.h"
#include <itkGradientDescentOptimizerv4.h>
#include <itkQuasiNewtonOptimizerv4.h>
#include <itkRegistrationParameterScalesFromIndexShift.h>
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
        void setupOptimizer(itk::GradientDescentOptimizerv4* opti);
        void setupOptimizer(itk::QuasiNewtonOptimizerv4* opti);
        void setupOptimizer(itk::FRPROptimizer* opti);

        void setupMetric(itk::MattesMutualInformationImageToImageMetricv4<RealImage, RealImage>* metric);
        void setupMetric(itk::EntropyImageToImageMetricv4<RealImage, RealImage>* metric);
        void setupMetric(itk::MeanSquaresImageToImageMetricv4<RealImage, RealImage>* metric);


        void setupTransform(itk::BSplineTransform<double,2,4>* transform);

    private:
        RealImage::Pointer _images[2];
        RealImage::Pointer _patches[2];
        RealImage::RegionType _patchRegion[2];
        QPolygonF _patchPolygon[2];

        typedef itk::TranslationTransform<double,2> TransformType;
//        typedef itk::Rigid2DTransform<double> TransformType;
//        typedef itk::BSplineTransform<double,2,4> TransformType;
        TransformType::Pointer _transform;

        typedef itk::GradientDescentOptimizerv4 OptimizerType;
        OptimizerType::Pointer _optimizer;

        typedef itk::MattesMutualInformationImageToImageMetricv4<RealImage, RealImage> CostFunctionType;
        CostFunctionType::Pointer _costFunc;

        typedef OptimizerType::ParametersType ParametersType;

        typedef itk::RegistrationParameterScalesFromIndexShift<CostFunctionType> ScaleEstimatorType;
        ScaleEstimatorType::Pointer _estimator;
    };
}
#endif /* defined(__ParticleGuidedRegistration__piPatchTracking__) */