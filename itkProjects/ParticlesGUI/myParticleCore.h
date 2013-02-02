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

namespace pi {

    class ParticleConstraint;
    class ImageContext;

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
        void AddForce(double* ff, double alpha = 1);
        void SubForce(double* ff, double alpha = 1);
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
        FieldTransformType::Pointer m_Transform;
        FieldTransformType::Pointer m_InverseTransform;

        ParticleSubject() : m_SubjId(-1) { }
        ParticleSubject(int subjid, int npoints);
        ~ParticleSubject();

        inline const int GetNumberOfPoints() const { return m_Particles.size(); }
        inline int GetNumberOfPoints() { return m_Particles.size(); }
        
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
    typedef boost::numeric::ublas::vector<ParticleSubject> ParticleSubjectArray;

    class InternalForce {
    public:
        InternalForce() : m_TimeStep(0.1) {}
        ~InternalForce() {}
        void ComputeForce(ParticleSubjectArray& shapes);
        void ComputeForce(Particle& a, Particle& b);
    private:
        double m_TimeStep;
    };

    class EnsembleForce {
    public:
        EnsembleForce();
        ~EnsembleForce();
        void SetImageContext(ImageContext* context);
        void ComputeEnsembleForce(ParticleSubjectArray& shapes);
        void ComputeImageForce(ParticleSubjectArray& shapes);
    private:
        ImageContext* m_ImageContext;
        ParticleSubject m_MeanShape;
        void ComputeMeanShape(ParticleSubjectArray& shapes);
    };

    class ImageContext {
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
    };


    class ParticleSystem {
    public:
        ParticleSystem(): m_NumParticlesPerSubject(300),
            m_ParticleConstraint(NULL) {
        }
        ~ParticleSystem() {
        }
        int GetNumberOfSubjects();
        int GetNumberOfParticles();
        ParticleSubjectArray& GetSubjects();
        ImageContext& GetImageContext();

        void LoadLabels(StringVector labels);
        void RunPreprocessing(std::string outputName);
        void LoadPreprocessing(std::string outputName);

        // run correspondence process
        void Run();
        void PrepareSystem(ParticleSubjectArray& shapes);
        void UpdateSystem(ParticleSubjectArray& shapes, double dt);


        void LoadSystem(std::string filename, int cmd = 0);
        void LoadSystem(std::string filename, ParticleSubjectArray& shapes, int cmd);
        void SaveSystem(std::string filename) { SaveSystem(filename, m_Subjects); }
        void SaveSystem(std::string filename, ParticleSubjectArray& shapes);

        inline ParticleSubject& operator[](int j) { return m_Subjects[j];
        }
        inline const ParticleSubject& operator[](int j) const {
            return m_Subjects[j];
        }
        
    private:
        int m_NumParticlesPerSubject;
        ImageContext m_ImageContext;
        ParticleSubjectArray m_Subjects;
        ParticleConstraint* m_ParticleConstraint;
        std::string m_TrackingOutputPattern;
    };
}

#endif
