//
//  myParticleForces.h
//  ParticlesGUI
//
//  Created by Joohwi Lee on 1/26/13.
//
//

#ifndef __ParticlesGUI__myParticleForces__
#define __ParticlesGUI__myParticleForces__

#include <iostream>
#include "piParticleCore.h"
#include "itkGradientImageFilter.h"

namespace pi {
    typedef itk::GradientImageFilter<DoubleImage> GradientFilterType;
    typedef GradientFilterType::OutputImageType GradientImage;
    typedef GradientFilterType::OutputPixelType GradientPixel;

    class InternalForce {
    public:
        InternalForce() {}
        ~InternalForce() {}
        void ComputeForce(ParticleSubject& subj);
        void ComputeForce(ParticleSubjectArray& shapes);
        void ComputeForce(Particle& a, Particle& b);
    };

    class EntropyInternalForce {
    public:
        EntropyInternalForce() {}
        ~EntropyInternalForce() {}

        void ComputeForce(ParticleSubject& subj);
        void ComputeForce(ParticleSubjectArray& subjs);
        void ComputeForce(Particle& a, Particle& b);
    };

    class EnsembleForce {
    public:
        EnsembleForce(double coeff);
        ~EnsembleForce();
        void SetImageContext(ImageContext* context);
        void ComputeEnsembleForce(ParticleSystem& system);
        void ComputeImageForce(ParticleSystem& system);
    private:
        ImageContext* m_ImageContext;
        ParticleSubject m_MeanShape;
        void ComputeMeanShape(ParticleSubjectArray& shapes);
        double m_Coeff;
    };


    template <int N>
    class ParticleAttribute3D {
    public:
        const static int NATTRS = N*N*N;
        double f[__Dim];
        double F[__Dim];
        double x[NATTRS];
        double y[NATTRS];
        double g[NATTRS][__Dim];
    };

    template <int N>
    class ParticleAttribute2D {
    public:
        const static int NATTRS = N*N;
        double f[__Dim];
        double F[__Dim];
        double x[NATTRS];
        double y[NATTRS];
        double g[NATTRS][__Dim];
    };

#ifdef DIMENSION3
    typedef ParticleAttribute3D<3> Attr;
#else
    typedef ParticleAttribute2D<3> Attr;
#endif

    class IntensityForce {
    public:
        IntensityForce(double coeff);
        ~IntensityForce();
        void SetImageContext(ImageContext* context);
        void SetWorkAtWarpedSpace(bool check);
        void ComputeIntensityForce(ParticleSystem* system);
    private:
        void ComputeAttributes(ParticleSystem* system);
        DoubleImageVector warpedImages;
        std::vector<GradientImage::Pointer> gradImages;
        typedef boost::numeric::ublas::matrix<Attr> AttrMatrix;
        AttrMatrix m_attrs;
        VNLMatrix m_attrsMean;
        ImageContext* m_ImageContext;
        double m_Coeff;
        bool m_WorkAtWarpedSpace;
    };
}

#endif /* defined(__ParticlesGUI__myParticleForces__) */
