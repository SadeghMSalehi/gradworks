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
#include "piEntropyComputer.h"


#ifndef ATTR_SIZE
#define ATTR_SIZE 3
#endif

class vtkPolyData;

namespace pi {

    class InternalForce {
    public:
        bool heteroForce;
        InternalForce(): heteroForce(false) {}
        ~InternalForce() {}
        void ComputeForce(ParticleSubject& subj);
        void ComputeForce(ParticleSubjectArray& shapes);
        void ComputeForce(Particle& a, Particle& b);
        DataReal repulsionSigma;
    };

    class EntropyInternalForce {
    public:
        EntropyInternalForce(): useMultiPhaseForce(false),repulsionSigma(5), repulsionCutoff(repulsionSigma*5),friendSigma(3), friendCutoff(friendSigma*5),  useAdaptiveSampling(false), maxKappa(3), coeff(1) {}
        ~EntropyInternalForce() {}

        void ComputeForce(ParticleSubject& subj);
        void ComputeForce(ParticleSubjectArray& subjs);
        void ComputeForce(Particle& a, Particle& b);

        bool useMultiPhaseForce;
        DataReal repulsionSigma;
        DataReal repulsionCutoff;
        DataReal friendSigma;
        DataReal friendCutoff;

        bool useAdaptiveSampling;
        DataReal maxKappa;
        DataReal coeff;

    private:
        void ComputeHomoForce(ParticleSubject& sub);
        void ComputeHeteroForce(ParticleSubject& sub);
    };

    

    class EnsembleForce {
    public:
        DataReal coeff;
        EnsembleForce();
        ~EnsembleForce();

        void SetImageContext(ImageContext* context);
        void ComputeEnsembleForce(ParticleSystem& system);
        void ComputeImageForce(ParticleSystem& system);
    private:
        ImageContext* m_ImageContext;

        bool useBSplineAlign;
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
    typedef ParticleAttribute3D<ATTR_SIZE> Attr;
#else
    typedef ParticleAttribute2D<ATTR_SIZE> Attr;
#endif

    class IntensityForce {
    public:
        typedef boost::numeric::ublas::matrix<Attr> AttrMatrix;

        DataReal coeff;
        bool useGaussianGradient;
        DataReal gaussianSigma;
        bool useAttributesAtWarpedSpace;

    public:
        IntensityForce();
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
    };
}

#endif /* defined(__ParticlesGUI__myParticleForces__) */
