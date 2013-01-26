//
//  ImageViewManager.h
//  laplacePDE
//
//  Created by Joohwi Lee on 11/13/12.
//
//

#ifndef __myParticlesCore__
#define __myParticlesCore__

#include "boost/numeric/ublas/vector.hpp"
#include "vnlCommon.h"
#include "itkImage.h"
#include "itkImageIO.h"
#include "itkOffset.h"
#include "vector"
#include "myImageDef.h"
#include "myImplicitSurfaceConstraint3D.h"

namespace my {

    class LabelContext {
    public:
        void LoadLabel(std::string filename);
        void ComputeIntersection();
        void ComputeDistanceMaps();
        LabelImage::Pointer GetLabel(int j);
        LabelImage::Pointer GetIntersection();
        OffsetImage GetDistanceMap(int j);
        ImplicitSurfaceConstraint* GetConstraint();
        
    private:
        LabelImage::Pointer m_Intersection;
        LabelVectors m_LabelImages;
        OffsetImageVectors m_DistanceMaps;
        ImplicitSurfaceConstraint m_Constraint;
    };

    const static double g_TimeStep = 0.1;

    class Particle {
    public:
        int subj;
        int idx;

        // the position x and the transformed point y
        double x[4];
        double y[4];

        // the current velocity v and the force f
        double v[4];
        double f[4];

        // the density and pressure of a particle
        double density;
        double pressure;

        Particle();
        ~Particle();

        inline void Set(LabelImage::IndexType idx) {
            x[0] = idx[0];
            x[1] = idx[1];
            x[2] = idx[2];
            x[3] = 1;
        }

        void Sub(const Particle& p, double* nx);
        void AddForce(double* ff);
        void SubForce(double* ff);
        double Dist2(const Particle& p);
        void UpdateVelocity(double *vv);
        void UpdateForce(double *ff);
        void UpdateSystem();

        Particle& operator=(const Particle& other);
    };


    typedef boost::numeric::ublas::vector<Particle> ParticleArray;

    class ParticleShape {
    public:
        int m_SubjId;
        int m_nPoints;
        ParticleArray m_Particles;
        
        ParticleShape(int subjid, int npoints);
        ~ParticleShape();

        void InitializeRandomPoints(LabelImage::Pointer labelImage);
        void Initialize(const ParticleArray& array);
        void ApplyMatrix(VNLMatrix& mat);
        void TransformX2Y(TransformType* transform);
        void TransformY2X(TransformType* transform);
        void UpdateInternalForce();
        void UpdateInternalConstraint();

        inline Particle& operator[](int i) {
            return m_Particles[i];
        }

        inline const Particle& operator[](int i) const {
            return m_Particles[i];
        }
    };
    typedef boost::numeric::ublas::vector<ParticleShape> ParticleShapeArray;

    class InternalForce {
    public:
        InternalForce() : m_TimeStep(0.1) {}
        ~InternalForce() {}
        void ComputeForce(Particle& a, Particle& b, double *ff);
    private:
        double m_TimeStep;
    };
    
    class EnsembleForce {
    public:
        EnsembleForce() {};
        ~EnsembleForce();
        
        void ComputeForce(ParticleShapeArray& shapes);
    };


    class ParticleSystem {
    public:
        void LoadShapes();
        void UpdateStep(double dt);
        void UpdateEnsembleForce();
    private:
        ParticleShapeArray m_Shapes;
    };
}

#endif