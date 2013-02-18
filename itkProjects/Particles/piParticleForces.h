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


#ifndef ATTR_DIM
#define ATTR_DIM 3
#endif

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
    typedef ParticleAttribute3D<ATTR_DIM> Attr;
#else
    typedef ParticleAttribute2D<ATTR_DIM> Attr;
#endif

    class IntensityForce {
    public:
        typedef boost::numeric::ublas::matrix<Attr> AttrMatrix;

        DataReal coefficient;
        bool useGaussianGradient;
        DataReal gaussianSigma;
        bool useAttributesAtWarpedSpace;

    public:
        IntensityForce(DataReal coeff);
        ~IntensityForce();
        
        void SetImageContext(ImageContext* context);
        void ComputeIntensityForce(ParticleSystem* system);

        // attrmatrix has S subject rows and N point columns
        void NormalizeAttributes(AttrMatrix& attrs);
        bool ComputeCovariance(AttrMatrix& attrs, int pointIdx, VNLDoubleMatrix& cov, double alpha = 1);
        void ComputeGradient(AttrMatrix& attrs, VNLDoubleMatrix& invC, int pointIdx, bool useDual);

    private:
        void ComputeAttributes(ParticleSystem* system);
        RealImageVector warpedImages;
        std::vector<GradientImage::Pointer> gradImages;
        AttrMatrix m_attrs;
        VNLMatrix m_attrsMean;
        ImageContext* m_ImageContext;
        DataReal m_Coeff;
    };
}

#endif /* defined(__ParticlesGUI__myParticleForces__) */
