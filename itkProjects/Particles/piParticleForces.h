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

    class InternalForce {
    public:
        InternalForce() {}
        ~InternalForce() {}
        void ComputeForce(ParticleSubject& subj);
        void ComputeForce(ParticleSubjectArray& shapes);
        void ComputeForce(Particle& a, Particle& b);
        DataReal repulsionSigma;
    };

    class EntropyInternalForce {
    public:
        EntropyInternalForce() :repulsionSigma(3), repulsionCutoff(repulsionSigma*5), useAdaptiveSampling(false), maxKappa(3) {}
        ~EntropyInternalForce() {}

        void ComputeForce(ParticleSubject& subj);
        void ComputeForce(ParticleSubjectArray& subjs);
        void ComputeForce(Particle& a, Particle& b);

        DataReal repulsionSigma;
        DataReal repulsionCutoff;
        bool useAdaptiveSampling;
        DataReal maxKappa;
    };

    class EnsembleForce {
    public:
        EnsembleForce(DataReal coeff);
        ~EnsembleForce();
        void SetImageContext(ImageContext* context);
        void ComputeEnsembleForce(ParticleSystem& system);
        void ComputeImageForce(ParticleSystem& system);
    private:
        ImageContext* m_ImageContext;
        ParticleSubject m_MeanShape;
        void ComputeMeanShape(ParticleSubjectArray& shapes);
        DataReal m_Coeff;
    };


    template <int N>
    class ParticleAttribute3D {
    public:
        const static int NATTRS = N*N*N;
        DataReal f[__Dim];
        DataReal F[__Dim];
        DataReal x[NATTRS];
        DataReal y[NATTRS];
        DataReal g[NATTRS][__Dim];
    };

    template <int N>
    class ParticleAttribute2D {
    public:
        const static int NATTRS = N*N;
        DataReal f[__Dim];
        DataReal F[__Dim];
        DataReal x[NATTRS];
        DataReal y[NATTRS];
        DataReal g[NATTRS][__Dim];
    };

#ifdef DIMENSION3
    typedef ParticleAttribute3D<3> Attr;
#else
    typedef ParticleAttribute2D<3> Attr;
#endif

    class IntensityForce {
    public:
        DataReal coefficient;
        bool useGaussianGradient;
        DataReal gaussianSigma;
        bool useAttributesAtWarpedSpace;

    public:
        IntensityForce(DataReal coeff);
        ~IntensityForce();
        
        void SetImageContext(ImageContext* context);
        void ComputeIntensityForce(ParticleSystem* system);
    private:
        void ComputeAttributes(ParticleSystem* system);
        DoubleImageVector warpedImages;
        std::vector<GradientImage::Pointer> gradImages;
        typedef boost::numeric::ublas::matrix<Attr> AttrMatrix;
        AttrMatrix m_attrs;
        VNLMatrix m_attrsMean;
        ImageContext* m_ImageContext;
        DataReal m_Coeff;
    };
}

#endif /* defined(__ParticlesGUI__myParticleForces__) */
