//
//  piParticleTools.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 2/16/13.
//
//

#ifndef __ParticleGuidedRegistration__piParticleTools__
#define __ParticleGuidedRegistration__piParticleTools__

#include <iostream>
#include "piImageDef.h"
#include "piParticleCore.h"
#include "piParticleBSpline.h"
#include "itkResampleImageFilter.h"
#include "itkCompositeTransform.h"

namespace pi {

    class ParticleTools {
    public:
        ParticleTools(Options& o, StringVector& a);
        void runCoverLabel();
        void particle2mat();
    private:
        Options& opts;
        StringVector& args;
    };

    template <class T>
    void MarkAtImage(T& data, int n, LabelImage::Pointer image, LabelPixel val) {
        for (int i = 0; i < n; i++) {
            Particle& pi = data[i];
            IntIndex idx;
            fordim (k) {
                idx[k] = pi.x[k] + 0.5;
            }
            image->SetPixel(idx, val);
        }
    }

    
    template <class S>
    typename S::Pointer warp_image(ParticleSubject& fixed, ParticleSubject& moving, typename S::Pointer m, typename S::Pointer ref, bool useNN, bool noRigidAlign, bool onlyRigidAlign) {

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
        cout << affineParams << endl;

        ParticleBSpline bsp;
        if (ref.IsNull()) {
            ref = m;
        }
        if (!onlyRigidAlign) {
            if (!noRigidAlign) {
                moving.TransformX2Y();
                bsp.EstimateTransform<ParticleYCaster,S>(fixed, moving, fixed.GetNumberOfPoints(), ref);
            } else {
                bsp.EstimateTransform<ParticleXCaster,S>(fixed, moving, fixed.GetNumberOfPoints(), ref);
            }
        }

        typedef itk::CompositeTransform<PointReal,__Dim> CompositeTransformType;
        typename CompositeTransformType::Pointer transform = CompositeTransformType::New();
        //        if (!noRigidAlign) {
        //            transform->AddTransform(affineTransform->GetInverseTransform());
        //        }
        if (!onlyRigidAlign) {
            transform->AddTransform(bsp.GetTransform());
        }
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
#endif /* defined(__ParticleGuidedRegistration__piParticleTools__) */
