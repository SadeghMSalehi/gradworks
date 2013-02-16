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
#include "boost/numeric/ublas/matrix.hpp"
#include "vector"
#include "piImageDef.h"
//#include "vnlCommon.h"
#include "itkImage.h"
#include "itkImageIO.h"
#include "itkOffset.h"
#include "piOptions.h"
#include "iostream"

namespace pi {

    class ParticleConstraint;
    class ImageContext;
    class ParticleSystem;

    class Particle {
    public:
        int subj;
        int idx;
        double t;

        // the position x and the transformed point y
        double x[4];
        double y[4];

        // the current velocity v and the force f
        double v[4];
        double f[4];

        // temporary status
        double w[4];

        // the density and pressure of a particle
        double density;
        double pressure;

        Particle();
        ~Particle();

        inline void Set(LabelImage::IndexType idx) {
            fordim(i) { x[i] = idx[i]; }
            x[__Dim] = 1;
        }
        
        void Zero();
        void Sub(const Particle& p, double* nx);
        void AddForce(const double* ff, double alpha = 1);
        void SubForce(const double* ff, double alpha = 1);
        double Dist2(const Particle& p);
        void UpdateVelocity(double *vv);
        void UpdateForce(double *ff);
        void UpdateSystem(double dt);

        Particle& operator=(const Particle& other);
    };
    // utility operator overloading
    ostream& operator<<(ostream& os, const Particle& par);
    istream& operator>>(istream& is, Particle& par);

    typedef boost::numeric::ublas::vector<Particle> ParticleArray;
    typedef std::vector<Particle> ParticleVector;
    
    class ParticleSubject {
    public:
        int m_SubjId;
        string m_Name;
        ParticleArray m_Particles;


        AffineTransformType::Pointer m_AffineTransform;
        FieldTransformType::Pointer m_DeformableTransform;
        FieldTransformType::Pointer m_InverseDeformableTransform;

        ParticleSubject() : m_SubjId(-1) { }
        ParticleSubject(int subjid, int npoints);
        ~ParticleSubject();

        inline const int GetNumberOfPoints() const { return m_Particles.size(); }
        inline int GetNumberOfPoints() { return m_Particles.size(); }

        void Clear();
        void Zero();
        void NewParticles(int nPoints);
        void InitializeRandomPoints(LabelImage::Pointer labelImage);
        void Initialize(int subj, std::string name, int nPoints);
        void Initialize(int subj, std::string name, const ParticleSubject& shape);
        void Initialize(const ParticleArray& array);
        void ApplyMatrix(VNLMatrix& mat);
        void TransformX2Y(TransformType* transform);
        void TransformY2X(TransformType* transform);
        void TransformX2X(TransformType* transform);
        void TransformY2Y(TransformType* transform);
        void ReadParticlePositions(std::istream& is, int nPoints);
        void WriteParticlePositions(std::ostream& os);
        void ReadParticles(std::istream& is, int nPoints);
        void WriteParticles(std::ostream& os);

        inline Particle& operator[](int i) {
            return m_Particles[i];
        }

        inline const Particle& operator[](int i) const {
            return m_Particles[i];
        }
    };

    ostream& operator<<(ostream& os, const Particle& par);

    typedef boost::numeric::ublas::vector<ParticleSubject> ParticleSubjectArray;


    class ImageContext {
        friend class ParticleSystem;
    public:
        void LoadLabel(std::string filename);
        void LoadDoubleImage(std::string filename);
        void ComputeIntersection();
        void ComputeDistanceMaps();
        LabelImage::Pointer GetLabel(int j);
        DoubleImage::Pointer GetDoubleImage(int j);
        LabelImage::Pointer GetIntersection();
        void SetIntersection(LabelImage::Pointer intersection);
        OffsetImage GetDistanceMap(int j);
        StringVector& GetDoubleImageFileNames();
        StringVector& GetFileNames();
        LabelVector& GetLabelVector();
        DoubleImageVector& GetDoubleImageVector();
        void Clear();

    private:
        StringVector m_DoubleImageFileNames;
        StringVector m_FileNames;
        LabelImage::Pointer m_Intersection;
        LabelVector m_LabelImages;
        DoubleImageVector m_Images;
        OffsetImageVector m_DistanceMaps;
        std::string m_IntersectionOutput;
    };

    class ParticleSystem {
    public:
        ParticleSystem();
        ~ParticleSystem() {
        }
        
        int GetNumberOfSubjects();
        int GetNumberOfParticles();
        void InitializeSystem(Options& options);
        
        ParticleSubject& GetInitialSubject();
        void ComputeMeanSubject();
        const ParticleSubject& GetMeanSubject() const;
        ParticleSubjectArray& GetSubjects();
        
        Options& GetSystemOptions() {
            return m_Options;
        }

        inline ParticleSubject& operator[](int j) {
            return m_Subjects[j];
        }

        inline const ParticleSubject& operator[](int j) const {
            return m_Subjects[j];
        }

    private:
        ParticleSubject m_MeanSubject;
        ParticleSubject m_InitialSubject;
        ParticleSubjectArray m_Subjects;
        
        Options m_Options;
    };

}

#endif
