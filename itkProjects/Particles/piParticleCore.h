//
//  ImageViewManager.h
//  laplacePDE
//
//  Created by Joohwi Lee on 11/13/12.
//
//

#ifndef __myParticlesCore__
#define __myParticlesCore__
#ifndef Q_MOC_RUN
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#endif
#include <vector>
#include <string>
#include <iostream>
#include <itkOffset.h>
#include <vnl/vnl_matrix_fixed.h>
#include "piImageDef.h"
#include "piOptions.h"
#include "vtkSmartPointer.h"
#include "vtkPoints.h"
#include "vtkTransform.h"
#include "vtkMatrix4x4.h"
#include "piParticle.h"

class vtkPolyData;

namespace pi {

    typedef vtkSmartPointer<vtkTransform> vtkTransformType;
    typedef vtkSmartPointer<vtkPoints> vtkPointsType;
    typedef vtkSmartPointer<vtkMatrix4x4> vtkMatrixType;

    class ImageContext;
    class ParticleSystem;

    // utility classes
    class ParticleXCaster {
    public:
        inline float castSource(const Particle& p, int i) const { return p.x[i]; }
        inline float castTarget(const Particle& p, int i) const { return p.x[i]; }
    };

    class ParticleYCaster {
    public:
        inline float castSource(const Particle& p, int i) const { return p.y[i]; }
        inline float castTarget(const Particle& p, int i) const { return p.y[i]; }
    };

    class ParticleZCaster {
    public:
        inline float castSource(const Particle& p, int i) const { return p.z[i]; }
        inline float castTarget(const Particle& p, int i) const { return p.z[i]; }
    };
    
    class ParticleYZCaster {
    public:
        inline float castSource(const Particle& p, int i) const { return p.y[i]; }
        inline float castTarget(const Particle& p, int i) const { return p.z[i]; }
    };

    ostream& operator<<(ostream& os, const vtkTransformType& par);
    istream& operator>>(istream& is, vtkTransformType& par);

    typedef boost::numeric::ublas::vector<Particle> ParticleArray;

    void createParticles(ParticleVector&, int subj, int n);

    class ParticleSubject {
    public:
        int m_SubjId;
        std::string m_Name;
        ParticleArray m_Particles;
        DataReal m_Spacing[DIMENSIONS];


    public:
        // Image related member variables
        std::string m_ImageFile;
        RealImageVector m_Images;
        RealImage::Pointer m_WarpedImage;

        std::string m_LabelFile;
        LabelImage::Pointer m_Label;

        OffsetImage::Pointer m_DistanceMap;
        
        RealImage::Pointer kappaImage;
        LinearImageInterpolatorType::Pointer kappaSampler;

        // Label sampler for multi-phase internal force
        LabelImage::Pointer friendImage;
        NNLabelInterpolatorType::Pointer friendSampler;


        RealImage::Pointer realImage;
        LinearImageInterpolatorType::Pointer realSampler;

        GradientImage::Pointer gradImage;
        GradientInterpolatorType::Pointer gradSampler;

        StringVector m_RealImageFileNames;
        StringVector m_FileNames;

    public:
        // traansform-related members
        vtkTransformType alignment;
        vtkTransformType inverseAlignment;
        vtkPointsType pointscopy;
//        AffineTransformType::Pointer m_AffineTransform;

    public:
        // deformable transforms
        FieldTransformType::Pointer m_DeformableTransform;
        FieldTransformType::Pointer m_InverseDeformableTransform;

    public:
        // public methods
        ParticleSubject();
        ParticleSubject(int subjid, int npoints);
        ~ParticleSubject();

        inline const int size() const { return m_Particles.size(); }
        inline const int GetNumberOfPoints() const { return m_Particles.size(); }
        inline int GetNumberOfPoints() { return m_Particles.size(); }

        void Clear();
        void Zero();
        void NewParticles(int nPoints);
        void SetSpacing(RealImage::SpacingType spacing);

        // initialize points
        void InitializeRandomPoints(LabelImage::Pointer labelImage);
        void Initialize(int subj, std::string name, int nPoints);
        void Initialize(int subj, std::string name, const ParticleSubject& shape);
        void Initialize(const ParticleArray& array);
        void SyncPointsCopy();

        // alignment
        void ComputeAlignment(ParticleSubject& subj, bool useSimilarity = false);
        void ComputeDensity();

        void AlignmentTransformX2Y();


        void SortByCorrespondence();
        void SortByCorrespondenceScore();

        // apply transform to particles
        void TransformX2Y(TransformType* transform = NULL);
        void TransformY2X(TransformType* transform = NULL);
        void TransformX2X(TransformType* transform);
        void TransformY2Y(TransformType* transform);
        void TransformY2Z(TransformType* transform = NULL);
        void TransformZ2Y(TransformType* transform = NULL);

        // conversion to index
        void ComputeIndexX(Particle& p, IntIndex& x);
        void ComputeIndexY(Particle& p, IntIndex& x);
        void ComputeIndexZ(Particle& p, IntIndex& x);
        void ComputeIndexX(Particle& p, RealIndex& x);
        void ComputeIndexY(Particle& p, RealIndex& x);
        void ComputeIndexZ(Particle& p, RealIndex& x);
        void SetFromIndex(IntIndex &x, Particle& p);

        void ContinuousIndexToPoint(ImagePoint& idx, ImagePoint& point);

        // particle IO
        void ReadParticlePositions(std::istream& is, int nPoints);
        void WriteParticlePositions(std::ostream& os);
        void ReadParticles(std::istream& is, int nPoints);
        void WriteParticles(std::ostream& os);
        void ReadAlignment(std::istream& is);
        void WriteAlignment(std::ostream& os);

        // image related functions
        void LoadLabel(std::string filename);
        void LoadImage(std::string filename);
        LabelImage::Pointer WarpLabelToMeanSpace();

        LabelImage::Pointer& GetLabel();
        void SetLabel(LabelImage::Pointer label);
        
        RealImage::Pointer& GetImage(int level = 0);
//
        // particle access
        inline Particle& operator[](int i) {
            return m_Particles[i];
        }

        inline const Particle& operator[](int i) const {
            return m_Particles[i];
        }


    };
    typedef boost::numeric::ublas::vector<ParticleSubject> ParticleSubjectArray;

 

    class ParticleSystem {
    public:
        DataReal currentTime;
        int currentIteration;

        // FIXME: image energy added
        VNLVector ImageEnergy;
        VNLVector NormalizedImageEnergy;

        ParticleSystem();
        ~ParticleSystem() {
        }

        const int size() const { return m_Subjects.size(); }
        const int points() const { return size() > 0 ? m_Subjects[0].m_Particles.size() : 0; }
        int GetNumberOfSubjects();
        int GetNumberOfParticles();

        void InitializeSystem(Options& options);
        void InitializeSystem(int nsubjs, int nparticles);
        void LoadKappaImages(Options& options);

        ParticleSubject& GetInitialSubject();
        ParticleSubject& InitializeMean();
        ParticleSubject& ComputeXMeanSubject();
        ParticleSubject& ComputeYMeanSubject();
        ParticleSubject& ComputeZMeanSubject();

        ParticleSubject& GetMeanSubject();
        ParticleSubjectArray& GetSubjects();

        void UseSingleParticle(int particleId);

        // initial intersection related function
        int ComputeIntersection();
        LabelImage::Pointer GetIntersection() {
            return m_Intersection;
        }

        int GetCurrentResolutionLevel();
        int GetMaximumResolutionLevel();
        void SetCurrentResolutionLevel(int level);

        void SetIntersection(LabelImage::Pointer intersection) {
            m_Intersection = intersection;
        }

        int SuggestNumberOfParticles();

        Options& GetSystemOptions() {
            return (*m_Options);
        }

        inline ParticleSubject& operator[](int j) {
            return m_Subjects[j];
        }

        inline const ParticleSubject& operator[](int j) const {
            return m_Subjects[j];
        }

        void RemoveParticle(int i);

        // public property to save intersection image to file
        std::string m_IntersectionOutput;

    private:
        ParticleSubject m_MeanSubject;
        ParticleSubject m_InitialSubject;
        ParticleSubjectArray m_Subjects;

        LabelImage::Pointer m_Intersection;

        int m_CurrentResolutionLevel;
        bool m_UseOriginalSpacing;

        Options* m_Options;
    };


    /*
    void export2vtk(ParticleSubject& sub, const char* vtkname, int field);
    vtkPolyData* convert2vtk(ParticleArray& parray);
     */
}

#endif
