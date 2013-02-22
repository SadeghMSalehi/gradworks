//
//  piImageWarp.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 2/21/13.
//
//

#include <stdio.h>

#include "piParticleCore.h"
#include "piParticleBSpline.h"

#include "itkResampleImageFilter.h"
#include "itkCompositeTransform.h"

namespace pi {
    template <class S>
    typename S::Pointer warp_image(ParticleSubject& fixed, ParticleSubject& moving, typename S::Pointer m, typename S::Pointer ref, bool useNN, bool noRigidAlign) {

        // We warp a given image m to m' based on a set of landmarks fixed and moving.
        // We do not assume that the landmarks are rigidly aligned.
        // Thus, we first estimate landmark-based rigid alignment T<fixed,moving> from fixed to moving
        // so that we can sample the fixed coordinate space from the moving space
        // After the rigid alignment, the bspline transform is estimated from the set of landmarks transformed to aligned space resulting in sampling from moving image

        fixed.ComputeAlignment(moving);
        fixed.AlignmentTransformX2Y();

        typedef itk::AffineTransform<PointReal, __Dim> AffineTransformType;
        AffineTransformType::Pointer affineTransform = AffineTransformType::New();
        AffineTransformType::ParametersType affineParams;
        affineParams.SetSize((__Dim+1)*(__Dim));
        for (int i = 0; i < __Dim; i++) {
            for (int j = 0; j < __Dim; j++) {
                affineParams[i*__Dim + j] = fixed.alignment->GetMatrix()->GetElement(i,j);
            }
            affineParams[__Dim*__Dim + i] = fixed.alignment->GetMatrix()->GetElement(i,3);
        }
        affineTransform->SetParameters(affineParams);
        
        ParticleBSpline bsp;
        if (ref.IsNull()) {
            ref = m;
        }
        bsp.EstimateTransform<ParticleYCaster,S>(fixed, moving, fixed.GetNumberOfPoints(), ref);

        typedef itk::CompositeTransform<PointReal,__Dim> CompositeTransformType;
        typename CompositeTransformType::Pointer transform = CompositeTransformType::New();
        transform->AddTransform(bsp.GetTransform());
        if (!noRigidAlign) {
            transform->AddTransform(affineTransform);
        }

        typedef itk::ResampleImageFilter<S,S> ResampleFilterType;
        typename ResampleFilterType::Pointer resampleFilter = ResampleFilterType::New();
        resampleFilter->SetInput(m);
        resampleFilter->SetTransform(transform);
        if (useNN) {
            typedef itk::NearestNeighborInterpolateImageFunction<S> InterpolatorType;
            typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
            resampleFilter->SetInterpolator(interpolator);
        }
        if (ref.IsNull()) {
            resampleFilter->SetReferenceImage(m);
        } else {
            resampleFilter->SetReferenceImage(ref);
        }
        resampleFilter->UseReferenceImageOn();
        resampleFilter->Update();
        return resampleFilter->GetOutput();
    }
}