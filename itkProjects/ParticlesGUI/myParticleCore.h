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

    class ParticleShape {
    public:
        int m_SubjId;
        int m_nPoints;
        string m_Name;
        ParticleArray m_Particles;
        FieldTransformType::Pointer m_Transform;
        
        ParticleShape() : m_SubjId(-1), m_nPoints(0) { Zero(); }
        ParticleShape(int subjid, int npoints);
        ~ParticleShape();

        void Zero();
        void NewParticles(int nPoints);
        void InitializeRandomPoints(LabelImage::Pointer labelImage);
        void Initialize(int subj, std::string name, const ParticleShape& shape);
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
    typedef boost::numeric::ublas::vector<ParticleShape> ParticleShapeArray;

    class InternalForce {
    public:
        InternalForce() : m_TimeStep(0.1) {}
        ~InternalForce() {}
        void ComputeForce(ParticleShapeArray& shapes);
        void ComputeForce(Particle& a, Particle& b);
    private:
        double m_TimeStep;
    };
    
    class EnsembleForce {
    public:
        EnsembleForce();
        ~EnsembleForce();
        void ComputeForce(ParticleShapeArray& shapes);
    private:
        ParticleShape m_MeanShape;
        void ComputeMeanShape(ParticleShapeArray& shapes);
    };

    class LabelContext {
    public:
        void LoadLabel(std::string filename);
        void ComputeIntersection();
        void ComputeDistanceMaps();
        LabelImage::Pointer GetLabel(int j);
        LabelImage::Pointer GetIntersection();
        OffsetImage GetDistanceMap(int j);
        StringVector& GetFileNames();
        void Clear();

    private:
        std::vector<std::string> m_FileNames;
        LabelImage::Pointer m_Intersection;
        LabelVector m_LabelImages;
        OffsetImageVector m_DistanceMaps;
    };

    
    class ParticleSystem {
    public:
        ParticleSystem(): m_ParticleConstraint(NULL) {
        }
        ~ParticleSystem() {
        }
        int GetNumberOfShapes();
        int GetNumberOfParticles();
        ParticleShapeArray& GetShapes();

        void LoadShapes(StringVector labels);
        void RunPreprocessing(std::string outputName);
        void LoadPreprocessing(std::string outputName);
        void UpdateStep(double dt);
        void PrepareSystem(ParticleShapeArray& shapes);
        void UpdateSystem(ParticleShapeArray& shapes, double dt);

				// run correspondence process
				void Run();

        void LoadStatus(std::string filename, int cmd = 0);
        void LoadStatus(std::string filename, ParticleShapeArray& shapes, int cmd);
        void SaveStatus(std::string filename) { SaveStatus(filename, m_Shapes); }
        void SaveStatus(std::string filename, ParticleShapeArray& shapes);

        inline ParticleShape& operator[](int j) { return m_Shapes[j];
        }
        inline const ParticleShape& operator[](int j) const {
            return m_Shapes[j];
        }

    private:
        LabelContext m_LabelContext;
        ParticleShapeArray m_Shapes;
        ParticleConstraint* m_ParticleConstraint;
    };
}

#endif
