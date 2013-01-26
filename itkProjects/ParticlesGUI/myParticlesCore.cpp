#include "myParticlesCore.h"

#include "myImplicitSurfaceConstraint3D.h"
#include "myParticleBSpline.h"

namespace my {

    // constructor
    // set every member variable as zero
    Particle::Particle() {
        for (int j = 0; j < 4; j++) {
            x[j] = y[j] = v[j] = f[j] = 0;
        }
        density = pressure = 0;
    }

    Particle::~Particle() {

    }

    void Particle::Sub(const Particle& p, double* f) {
        for3(i) {
            f[i] = x[i] - p.x[i];
        }
    }

    void Particle::AddForce(double* ff) {
        for3(i) {
            f[i] += ff[i];
        }
    }

    void Particle::SubForce(double* ff) {
        for3(i) {
            f[i] -= ff[i];
        }
    }

    double Particle::Dist2(const Particle& p) {
        double y[4];
        for4(i) {
            y[i] = x[i] - p.x[i];
        }
        return y[0]*y[0]+y[1]*y[1]+y[2]*y[2];
    }

    void Particle::UpdateForce(double *ff) {
        for4(i) {
            f[i] = ff[i];
        }
    }

    void Particle::UpdateVelocity(double *vv) {
        for4(i) {
            v[i] = vv[i];
        }
    }

    void Particle::UpdateSystem() {
        for3(i) {
            x[i] += g_TimeStep * v[i];
            v[i] += g_TimeStep * f[i];
        }
    }

    Particle& Particle::operator=(const Particle& other) {
        for4(i) {
            x[i] = other.x[i];
            y[i] = other.y[i];
            v[i] = other.v[i];
            f[i] = other.f[i];
        }
        return (*this);
    }

    ParticleShape::ParticleShape(int subjid, int npoints) : m_SubjId(subjid), m_nPoints(npoints) {

    }

    ParticleShape::~ParticleShape() {

    }

    void ParticleShape::InitializeRandomPoints(LabelImage::Pointer intersection) {
        // resize particles array
        m_Particles.resize(m_nPoints);

        // compute intersection by looping over region
        std::vector<LabelImage::IndexType> indexes;
        LabelImageIteratorType iter(intersection, intersection->GetBufferedRegion());
        for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
            LabelImage::IndexType idx = iter.GetIndex();
            LabelImage::PixelType pixel = iter.Value();
            if (pixel > 0) {
                indexes.push_back(idx);
            }
        }
        if (indexes.size() > 0) {
            std::random_shuffle(indexes.begin(), indexes.end());
        }

        
        for (int i = 0; i < m_nPoints; i++) {
            LabelImage::IndexType idx = indexes[i];
            m_Particles[i].Set(idx);
        }
    }


    void ParticleShape::Initialize(const ParticleArray& array) {
        if (m_nPoints != array.size()) {
            m_nPoints = array.size();
        }
        m_Particles.resize(m_nPoints);
        for (int i = 0; i < m_nPoints; i++) {
            
        }
    }

    void ParticleShape::ApplyMatrix(VNLMatrix &mat) {

    }

    void ParticleShape::ApplyTransform(TransformType* transform) {
        for (int i = 0; i < m_nPoints; i++) {
            TransformType::InputPointType inputPoint;
            for3 (j) {
                inputPoint[j] = m_Particles[i].x[j];
            }
            TransformType::OutputPointType outputPoint = transform->TransformPoint(inputPoint);
            for3 (j) {
                m_Particles[i].y[j] = outputPoint[j];
            }
        }
    }

    void ParticleShape::UpdateInternalForce() {
        InternalForce forceAlg;
        for (int i = 0; i < m_nPoints; i++) {
            for (int j = i+1; j < m_nPoints; j++) {
                double f[3];
                forceAlg.ComputeForce(m_Particles[i], m_Particles[j], f);
            }
        }
    }

    void ParticleSystem::LoadShapes() {

    }

    void ParticleSystem::UpdateStep(double dt) {
        for (int i = 0; i < m_Shapes.size(); i++) {
            m_Shapes[i].UpdateInternalForce();
            if (i > 0) {
                ParticleBSpline transform;
                transform.EstimateTransform(m_Shapes[0], m_Shapes[i]);
                m_Shapes[i].ApplyTransform(transform.GetTransform().GetPointer());
            }
        }
    }


    void InternalForce::ComputeForce(my::Particle &a, my::Particle &b, double* f) {
        double dist = std::sqrt(a.Dist2(b));
        const double sigma = 7 * 5;
        const double coeff = M_PI_2 / sigma;

        double dx[4];
        a.Sub(b, dx);

        if (dist > sigma) {
            return;
        }

        double rij = dist * coeff;
        if (rij > 0) {
            double sin2rij = std::sin(rij);
            sin2rij *= sin2rij;
            for3(i) {
                f[i] = dx[i] * (coeff * (1 - 1 / sin2rij));
            }
        }
        b.AddForce(f);
        a.SubForce(f);
    }

    void LabelContext::LoadLabel(std::string filename) {
        itkcmds::itkImageIO<LabelImage> io;
        LabelImage::Pointer image = io.ReadImageT(filename.c_str());
        m_LabelImages.push_back(image);
    }

    LabelImage::Pointer LabelContext::GetLabel(int j) {
        return m_LabelImages[j];
    }

    ImplicitSurfaceConstraint* LabelContext::GetConstraint() {
        return &m_Constraint;
    }

    void LabelContext::ComputeIntersection() {
        itkcmds::itkImageIO<LabelImage> io;
        LabelImage::Pointer intersection = io.NewImageT(m_LabelImages[0]);
        LabelImage::RegionType region = intersection->GetBufferedRegion();

        // set as member variable to reuse
        m_Intersection = intersection;

        // compute intersection by looping over region
        std::vector<LabelImage::IndexType> indexes;
        LabelImageIteratorType iter(intersection, region);
        for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
            LabelImage::IndexType idx = iter.GetIndex();
            LabelImage::PixelType pixel = 255;
            const int pixelThreshold = 255;
            for (int i = 0; i < m_LabelImages.size(); i++) {
                if (m_LabelImages[i]->GetPixel(idx) < pixelThreshold) {
                    pixel = 0;
                    break;
                }
            }
            if (pixel > 0) {
                intersection->SetPixel(idx, 255);
            }
        }
        io.WriteImageT("/tmpfs/intersection.nrrd", intersection);
    }

    void LabelContext::ComputeDistanceMaps() {
        m_Constraint.Clear();
        m_Constraint.SetImageList(m_LabelImages);
    }
}