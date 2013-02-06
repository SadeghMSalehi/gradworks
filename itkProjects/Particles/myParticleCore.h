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
//#include "vnlCommon.h"
#include "itkImage.h"
#include "itkImageIO.h"
#include "itkOffset.h"
#include "vector"
#include "myImageDef.h"

namespace pi {

    class ParticleConstraint;
    class ImageContext;
    class ParticleSystem;

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


    typedef boost::numeric::ublas::vector<Particle> ParticleArray;

    class ParticleSubject {
    public:
        int m_SubjId;
        string m_Name;
        ParticleArray m_Particles;
        CompositeTransformType::Pointer m_Transform;
        CompositeTransformType::Pointer m_InverseTransform;

        ParticleSubject() : m_SubjId(-1) { }
        ParticleSubject(int subjid, int npoints);
        ~ParticleSubject();

        inline const int GetNumberOfPoints() const { return m_Particles.size(); }
        inline int GetNumberOfPoints() { return m_Particles.size(); }

        void Clear();
        void Zero();
        void NewParticles(int nPoints);
        void InitializeRandomPoints(LabelImage::Pointer labelImage);
        void Initialize(int subj, std::string name, const ParticleSubject& shape);
        void Initialize(const ParticleArray& array);
        void ApplyMatrix(VNLMatrix& mat);
        void TransformX2Y(TransformType* transform);
        void TransformY2X(TransformType* transform);
        void UpdateSystem(double dt);

        inline Particle& operator[](int i) {
            return m_Particles[i];
        }

        inline const Particle& operator[](int i) const {
            return m_Particles[i];
        }
    };

    ostream& operator<<(ostream& os, const Particle& par);

    typedef boost::numeric::ublas::vector<ParticleSubject> ParticleSubjectArray;

    class InternalForce {
    public:
        InternalForce() {}
        ~InternalForce() {}
        void ComputeForce(ParticleSubjectArray& shapes);
        void ComputeForce(Particle& a, Particle& b);
    };

    class EnsembleForce {
    public:
        EnsembleForce(double coeff);
        ~EnsembleForce();
        void SetImageContext(ImageContext* context);
        void ComputeEnsembleForce(ParticleSubjectArray& shapes);
        void ComputeImageForce(ParticleSubjectArray& shapes);
    private:
        ImageContext* m_ImageContext;
        ParticleSubject m_MeanShape;
        void ComputeMeanShape(ParticleSubjectArray& shapes);
        double m_Coeff;
    };

    class IntensityForce {
    public:
        IntensityForce(double coeff);
        ~IntensityForce();
        void SetImageContext(ImageContext* context);
        void ComputeIntensityForce(ParticleSystem* system);
    private:
        ImageContext* m_ImageContext;
        double m_Coeff;
    };

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

    class ParticleSlice {
    public:
        typedef std::vector<Particle*> ParticlePointerVector;

        ParticleSlice(int dim, int nSubj, int nPoints) {
            m_SliceDim = dim;
            m_NumSubjects = nSubj;
            m_NumPoints = nPoints;
        }

        const ParticlePointerVector& Get(int slice, int subj);
        void Update(ParticleSubjectArray& shapes, LabelImage::Pointer refImage);

    private:
        int m_SliceDim;
        int m_NumSubjects;
        int m_NumPoints;
        
        boost::numeric::ublas::matrix<ParticlePointerVector> m_ParticlePointerMatrix;
    };

    class ParticleSystem {
    public:
        ParticleSystem();
        ~ParticleSystem() {
        }
        int GetNumberOfSubjects();
        int GetNumberOfParticles();

        void ComputeMeanSubject();
        const ParticleSubject& GetMeanSubject() const;
        ParticleSubjectArray& GetSubjects();
        ImageContext& GetImageContext();

        void LoadLabels(StringVector labels);
        void RunPreprocessing();

        // run correspondence process
        void Run();
        void UpdateSystem(ParticleSubjectArray& shapes, double dt);
        //        void PrepareSystem(ParticleSubjectArray& shapes);



        bool LoadSystem(std::string filename);
        void SaveSystem(std::string filename);
        
        inline ParticleSubject& operator[](int j) {
            return m_Subjects[j];
        }

        inline const ParticleSubject& operator[](int j) const {
            return m_Subjects[j];
        }
        
    private:
        double m_times[3];
        double m_timesPreprocessing[3];

        // flags for force selection
        bool m_InternalForceFlag;
        bool m_EnsembleForceFlag;
        bool m_IntensityForceFlag;
        
        double m_EnsembleCoeff;
        double m_IntensityCoeff;

        int m_NumParticlesPerSubject;
        ImageContext m_ImageContext;
        ParticleSubject m_MeanSubject;
        ParticleSubjectArray m_Initial;
        ParticleSubjectArray m_Subjects;
        std::string m_TrackingOutputPattern;
    };
}

#endif
