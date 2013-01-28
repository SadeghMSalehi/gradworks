#include "myParticleCore.h"

#include "myParticleConstraint.h"
#include "myParticleBSpline.h"
#include "iostream"
#include "sstream"

namespace pi {

    // constructor
    // set every member variable as zero
    Particle::Particle() {
        Zero();
    }

    Particle::~Particle() {

    }

    void Particle::Zero() {
        for4(j) {
            x[j] = y[j] = v[j] = f[j] = 0;
        }
        density = pressure = 0;
    }

    void Particle::Sub(const Particle& p, double* d) {
        fordim(i) {
            d[i] = x[i] - p.x[i];
        }
    }

    void Particle::AddForce(double* ff, double alpha) {
        fordim(i) {
            f[i] += (alpha * ff[i]);
        }
    }

    void Particle::SubForce(double* ff, double alpha) {
        fordim(i) {
            f[i] -= (alpha * ff[i]);
        }
    }

    double Particle::Dist2(const Particle& p) {
        double d[__Dim];
        fordim(i) {
            d[i] = x[i] - p.x[i];
        }
        double dist2 = 0;
        fordim(k) {
            dist2 += (d[k]*d[k]);
        }
        return dist2;
    }

    void Particle::UpdateForce(double *ff) {
        fordim(i) {
            f[i] = ff[i];
        }
    }

    void Particle::UpdateVelocity(double *vv) {
        fordim(i) {
            v[i] = vv[i];
        }
    }

    void Particle::UpdateSystem(double dt) {
        fordim(i) {
            x[i] += dt * v[i];
            v[i] += dt * f[i];
        }
    }

    Particle& Particle::operator=(const Particle& other) {
        for4(i) {
            x[i] = other.x[i];
            y[i] = other.y[i];
            v[i] = other.v[i];
            f[i] = other.f[i];
            density = other.density;
            pressure = other.pressure;
        }
        return (*this);
    }

    ParticleShape::ParticleShape(int subjid, int npoints) : m_SubjId(subjid), m_nPoints(npoints) {

    }

    ParticleShape::~ParticleShape() {

    }

    void ParticleShape::Zero() {
        for (int i = 0; i < m_nPoints; i++) {
            m_Particles[i].Zero();
        }
    }

    void ParticleShape::NewParticles(int n) {
        m_nPoints = n;
        m_Particles.resize(n);
        Zero();
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

    void ParticleShape::Initialize(int subj, std::string name, const ParticleShape& shape) {
        m_SubjId = subj;
        m_Name = name;
        m_nPoints = shape.m_Particles.size();
        m_Particles.resize(m_nPoints);
        for (int i = 0; i < m_nPoints; i++) {
            m_Particles[i] = shape.m_Particles[i];
        }
    }

    void ParticleShape::Initialize(const ParticleArray& array) {
        if (m_nPoints != array.size()) {
            m_nPoints = array.size();
        }
        m_Particles.resize(m_nPoints);
        for (int i = 0; i < m_nPoints; i++) {
            m_Particles[i] = array[i];
        }
    }

    void ParticleShape::ApplyMatrix(VNLMatrix &mat) {

    }

    void ParticleShape::TransformX2Y(TransformType* transform) {
        if (transform != NULL) {
            for (int i = 0; i < m_nPoints; i++) {
                TransformType::InputPointType inputPoint;
                fordim (j) {
                    inputPoint[j] = m_Particles[i].x[j];
                }
                TransformType::OutputPointType outputPoint = transform->TransformPoint(inputPoint);
                fordim (j) {
                    m_Particles[i].y[j] = outputPoint[j];
                }
            }
        } else {
            for (int i = 0; i < m_nPoints; i++) {
                for4(j) {
                    m_Particles[i].y[j] = m_Particles[i].x[j];
                }
            }
        }
    }
    
    void ParticleShape::TransformY2X(TransformType* transform) {
        if (transform != NULL) {
            for (int i = 0; i < m_nPoints; i++) {
                TransformType::InputPointType inputPoint;
                fordim (j) {
                    inputPoint[j] = m_Particles[i].y[j];
                }
                TransformType::OutputPointType outputPoint = transform->TransformPoint(inputPoint);
                fordim (j) {
                    m_Particles[i].x[j] = outputPoint[j];
                }
            }
        } else {
            for (int i = 0; i < m_nPoints; i++) {
                for4(j) {
                    m_Particles[i].x[j] = m_Particles[i].y[j];
                }
            }
        }
    }

    void ParticleShape::UpdateSystem(double dt) {
        for (int i = 0; i < m_nPoints; i++) {
            m_Particles[i].UpdateSystem(dt);
        }
    }

    int ParticleSystem::GetNumberOfShapes() {
        return m_Shapes.size();
    }

    int ParticleSystem::GetNumberOfParticles() {
        if (m_Shapes.size() > 0) {
            return m_Shapes[0].m_nPoints;
        }
        return 0;
    }

    ParticleShapeArray& ParticleSystem::GetShapes() {
        return m_Shapes;
    }

    void ParticleSystem::LoadShapes(StringVector files) {
        for (int i = 0; i < files.size(); i++) {
            m_LabelContext.LoadLabel(files[i]);
        }
    }

    void ParticleSystem::LoadPreprocessing(std::string filename) {
        ParticleShapeArray shapes;
        LoadStatus(filename, shapes, 1);
        StringVector& labelFiles = m_LabelContext.GetFileNames();
        const int nShapes = labelFiles.size();
        m_Shapes.resize(nShapes);
        for (int i = 0; i < nShapes; i++) {
            m_Shapes[i].Initialize(i, labelFiles[i], shapes[0]);
        }
    }

    void ParticleSystem::RunPreprocessing(std::string outputName) {
        double t0 = 0;
        double t1 = 10;
        double dt = 0.05;

        m_LabelContext.ComputeIntersection();
        
        ParticleShapeArray shapes(1);
        shapes[0].m_Name = "Intersection";
        shapes[0].m_nPoints = 300;
        shapes[0].InitializeRandomPoints(m_LabelContext.GetIntersection());

        ParticleConstraint initialConstraint;
        LabelVector labels;
        labels.push_back(m_LabelContext.GetIntersection());
        initialConstraint.SetImageList(labels);

        int k = 0;
        char trackName[128];

        for (double t = t0; t < t1; t += dt) {
            cout << "Time: " << t << endl;
            InternalForce internalForce;
            internalForce.ComputeForce(shapes);
            initialConstraint.ApplyConstraint(shapes);
            UpdateSystem(shapes, dt);
            sprintf(trackName, "preprocessing_%04d.txt", ++k);
            SaveStatus(trackName, shapes);
        }
        SaveStatus(outputName);
    }

    void ParticleSystem::UpdateStep(double dt) {
        InternalForce internalForce;
        EnsembleForce ensembleForce;
        internalForce.ComputeForce(m_Shapes);
        ensembleForce.ComputeForce(m_Shapes);
        m_ParticleConstraint->ApplyConstraint(m_Shapes);
    }

    void ParticleSystem::PrepareSystem(ParticleShapeArray& shapes) {
        if (shapes.size() < 1) {
            return;
        }
        const int nShapes = shapes.size();
        const int nPoints = shapes[0].m_Particles.size();
        for (int i = 0; i < nShapes; i++) {
            for (int j = 0; j < nPoints; j++) {
                fordim(k) {
                    shapes[i][j].f[k] = 0;
                }
            }
        }
    }

    void ParticleSystem::UpdateSystem(ParticleShapeArray& shapes, double dt) {
        const int nShapes = shapes.size();
        const int nPoints = shapes[0].m_Particles.size();
        for (int i = 0; i < nShapes; i++) {
            ParticleShape& iShape = shapes[i];
            for (int j = 0; j < nPoints; j++) {
                Particle& p = iShape[j];
                fordim(k) {
                    p.x[k] += (dt * p.v[k]);
                    p.v[k] += (dt * p.f[k]);
                    p.f[k] = 0;
                }
            }
        }
    }

    void ParticleSystem::LoadStatus(std::string filename, int cmd) {
        LoadStatus(filename, m_Shapes, cmd);
    }

    void ParticleSystem::LoadStatus(std::string filename, ParticleShapeArray& shapes, int cmd) {
        using namespace std;
        int nShapes = 0;
        int nParticles = 0;

        m_LabelContext.Clear();

        char cbuf[256];
        ifstream in(filename.c_str());

        string name;
        int value;

        in.getline(cbuf, sizeof(cbuf));
        {
            stringstream ss(cbuf);
            ss >> name >> value;
            for (int i = 0; i < value; i++) {
                in.getline(cbuf, sizeof(cbuf));
                if (cmd != 0) {
                    m_LabelContext.LoadLabel(cbuf);
                }
            }
        }

        in.getline(cbuf, sizeof(cbuf));
        {
            shapes.clear();
            stringstream ss(cbuf);
            ss >> name >> nShapes;
            shapes.resize(nShapes);
            for (int i = 0; i < nShapes; i++) {
                in.getline(cbuf, sizeof(cbuf));
                shapes[i].m_Name = cbuf;
            }
        }

        in.getline(cbuf, sizeof(cbuf));
        {
            stringstream ss(cbuf);
            ss >> name >> nParticles;

            cout << "Shapes: " << nShapes << ", Particles = " << nParticles << endl;
            for (int i = 0; i < nShapes; i++) {
                shapes[i].NewParticles(nParticles);
                for (int j = 0; j < nParticles; j++) {
                    Particle& p = shapes[i][j];
                    in.getline(cbuf, sizeof(cbuf));
                    {
                        stringstream ss(cbuf);
                        for4(k) { ss >> p.x[k]; }
                        for4(k) { ss >> p.y[k]; }
                        for4(k) { ss >> p.v[k]; }
                        for4(k) { ss >> p.f[k]; }
                        ss >> p.density;
                        ss >> p.pressure;
                    }
                }
            }
        }
    }

    void ParticleSystem::SaveStatus(std::string filename, ParticleShapeArray& shapes) {
        using namespace std;
        ofstream out(filename.c_str());
        StringVector& labels = m_LabelContext.GetFileNames();
        out << "LabelImages: " << labels.size() << endl;
        for (int i = 0; i < labels.size(); i++) {
            out << labels[i] << endl;
        }

        const int nShapes = shapes.size();
        const int nParticles = shapes[0].m_Particles.size();
        out << "Shapes: " << shapes.size() << endl;
        for (int i = 0; i < nShapes; i++) {
            out << m_Shapes[i].m_Name << endl;
        }
        out << "Particles: " << shapes[0].m_Particles.size() << endl;
        for (int i = 0; i < nShapes; i++) {
            for (int j = 0; j < nParticles; j++) {
                Particle& p = shapes[i][j];
                for4(k) { out << p.x[k] << " "; }
                for4(k) { out << p.y[k] << " "; }
                for4(k) { out << p.v[k] << " "; }
                for4(k) { out << p.f[k] << " "; }
                out << p.density << " ";
                out << p.pressure << " ";
                out << endl;
            }
        }
        out.close();
    }

    void LabelContext::Clear() {
        m_FileNames.clear();
        m_LabelImages.clear();
        m_DistanceMaps.clear();
    }
    
    void LabelContext::LoadLabel(std::string filename) {
        itkcmds::itkImageIO<LabelImage> io;
        LabelImage::Pointer image = io.ReadImageT(filename.c_str());
        m_LabelImages.push_back(image);
        m_FileNames.push_back(filename);
    }

    LabelImage::Pointer LabelContext::GetLabel(int j) {
        return m_LabelImages[j];
    }

    LabelImage::Pointer LabelContext::GetIntersection() {
        return m_Intersection;
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

    StringVector& LabelContext::GetFileNames() {
        return m_FileNames;
    }
}
