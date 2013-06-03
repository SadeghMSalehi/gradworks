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
        EntropyInternalForce(): useMultiPhaseForce(false),repulsionSigma(5), repulsionCutoff(repulsionSigma*5),friendSigma(0.45), friendCutoff(friendSigma*5),  useAdaptiveSampling(false), maxKappa(3*0.15), coeff(1) {}
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


    class NeighborSampler {
    public:
        typedef RealImage::PixelType PixelType;
        typedef RealImage::IndexType IndexType;
        typedef RealImage::PointType PointType;
        typedef RealImage::RegionType RegionType;
        typedef RealImage::SizeType SizeType;

        NeighborSampler(RegionType region, RealImage::Pointer image);
        ~NeighborSampler();

        void setSampleRegion(RegionType& region, RealImage* reference);
        
        void sampleValues(LinearImageInterpolatorType* interp, Particle& particle, ParticleAttribute &attr);
        void sampleGradients(GradientInterpolatorType* interp, Particle& particle, ParticleAttribute &attr);

        void createSampleIndexes2(IndexType& startIdx);
        void createSampleIndexes3(IndexType& startIdx);
        
    private:
        RealImage::SizeType _regionSize;

        int _numberOfSamples;
        std::vector<IndexType> _indexes;
        std::vector<PointType> _points;
    };

    class IntensityForce {
    public:
        typedef boost::numeric::ublas::matrix<ParticleAttribute> AttrMatrix;

        DataReal coeff;
        bool useGaussianGradient;
        DataReal gaussianSigma;
        bool useAttributesAtWarpedSpace;
        bool computeForwardTransform;

    public:
        IntensityForce();
        ~IntensityForce();
        
        void SetImageContext(ImageContext* context);
        void ComputeIntensityForce(ParticleSystem* system);

        void ComputeTransform(ParticleSystem* system);
        void SampleAttributes(ParticleSystem* system);

        
        // attrmatrix has S subject rows and N point columns
        void NormalizeAttributes(AttrMatrix& attrs);
        bool ComputeCovariance(AttrMatrix& attrs, int pointIdx, VNLDoubleMatrix& cov, double alpha = 1);
        void ComputeGradient(AttrMatrix& attrs, VNLDoubleMatrix& invC, int pointIdx, bool useDual);

        // access to particle attribute
        ParticleAttribute* GetAttribute(int imageId, int particleId);

        // deprecated function
        void ComputeAttributes(ParticleSystem* system);

    private:
        RealImageVector warpedImages;
        std::vector<GradientImage::Pointer> gradImages;
        AttrMatrix m_attrs;
        VNLMatrix m_attrsMean;
        ImageContext* m_ImageContext;
    };
}

#endif /* defined(__ParticlesGUI__myParticleForces__) */
