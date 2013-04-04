//
//  piFitCurve.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/3/13.
//
//

#ifndef __ParticleGuidedRegistration__piFitCurve__
#define __ParticleGuidedRegistration__piFitCurve__

#include <iostream>
#include "piImageDef.h"
#include "piParticleCore.h"
#include <itkBSplineScatteredDataPointSetToImageFilter.h>
#include <itkVectorLinearInterpolateImageFunction.h>

namespace pi {
    class CurveFitting {
    public:
        // definition for displacement field
        typedef itk::Vector<DataReal,1> ScalarPoint;
        typedef itk::PointSet<ScalarPoint,1> SparseScalarSet;
        typedef itk::Image<ScalarPoint,1> DenseScalarImage;
        typedef itk::BSplineScatteredDataPointSetToImageFilter<SparseScalarSet, DenseScalarImage> BSplineFitter;
        typedef itk::VectorLinearInterpolateImageFunction<DenseScalarImage> BSplineInterpolator;

        CurveFitting();
        ~CurveFitting();

        pi::ParticleVector& GetResult();
        void FitCurve(pi::ParticleVector& particles, bool closed = true);
        void ApproximateCurve(pi::ParticleVector& particles, bool closed = true);

    private:
        DenseScalarImage::Pointer OneDimensionalBsplineFitting(pi::ParticleVector& particles, int dim, bool closed = true);
        void ComputeChordLengthParameterization(pi::ParticleVector& particles);

    private:
        pi::ParticleVector _result;
    };
}
#endif /* defined(__ParticleGuidedRegistration__piFitCurve__) */
