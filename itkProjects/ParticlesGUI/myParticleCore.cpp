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

    ParticleSubject::ParticleSubject(int subjid, int npoints) : m_SubjId(subjid) {
        NewParticles(npoints);
    }

    ParticleSubject::~ParticleSubject() {

    }

    void ParticleSubject::Zero() {
        const int nPoints = GetNumberOfPoints();
        for (int i = 0; i < nPoints; i++) {
            m_Particles[i].Zero();
        }
    }

    void ParticleSubject::NewParticles(int n) {
        m_Particles.resize(n);
        Zero();
    }

    void ParticleSubject::InitializeRandomPoints(LabelImage::Pointer intersection) {
        // resize particles array
        const int nPoints = GetNumberOfPoints();

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


        for (int i = 0; i < nPoints; i++) {
            LabelImage::IndexType idx = indexes[i];
            m_Particles[i].Set(idx);
        }
    }

    void ParticleSubject::Initialize(int subj, std::string name, const ParticleSubject& shape) {
        m_SubjId = subj;
        m_Name = name;
        const int nPoints = shape.GetNumberOfPoints();
        m_Particles.resize(nPoints);
        for (int i = 0; i < nPoints; i++) {
            m_Particles[i] = shape.m_Particles[i];
        }
    }

    void ParticleSubject::Initialize(const ParticleArray& array) {
        const int nPoints = GetNumberOfPoints();
        if (nPoints != array.size()) {
            m_Particles.resize(nPoints);
        }

        for (int i = 0; i < nPoints; i++) {
            m_Particles[i] = array[i];
        }
    }

    void ParticleSubject::ApplyMatrix(VNLMatrix &mat) {

    }

    void ParticleSubject::TransformX2Y(TransformType* transform) {
        const int nPoints = GetNumberOfPoints();

        if (transform != NULL) {
            for (int i = 0; i < nPoints; i++) {
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
            for (int i = 0; i < nPoints; i++) {
                for4(j) {
                    m_Particles[i].y[j] = m_Particles[i].x[j];
                }
            }
        }
    }

    void ParticleSubject::TransformY2X(TransformType* transform) {
        const int nPoints = GetNumberOfPoints();

        if (transform != NULL) {
            for (int i = 0; i < nPoints; i++) {
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
            for (int i = 0; i < nPoints; i++) {
                for4(j) {
                    m_Particles[i].x[j] = m_Particles[i].y[j];
                }
            }
        }
    }

    void ParticleSubject::UpdateSystem(double dt) {
        const int nPoints = GetNumberOfPoints();

        for (int i = 0; i < nPoints; i++) {
            m_Particles[i].UpdateSystem(dt);
        }
    }

    int ParticleSystem::GetNumberOfSubjects() {
        return m_Subjects.size();
    }

    int ParticleSystem::GetNumberOfParticles() {
        if (m_Subjects.size() > 0) {
            return m_Subjects[0].GetNumberOfPoints();
        }
        return 0;
    }

    ParticleSubjectArray& ParticleSystem::GetSubjects() {
        return m_Subjects;
    }

    ImageContext& ParticleSystem::GetImageContext() {
        return m_ImageContext;
    }

    void ParticleSystem::LoadLabels(StringVector files) {
        for (int i = 0; i < files.size(); i++) {
            m_ImageContext.LoadLabel(files[i]);
        }
    }

    void ParticleSystem::LoadPreprocessing(std::string filename) {
        ParticleSubjectArray subjects;
        LoadSystem(filename, subjects, 1);
        StringVector& labelFiles = m_ImageContext.GetFileNames();
        const int nSubjects = labelFiles.size();
        m_Subjects.resize(nSubjects);
        for (int i = 0; i < nSubjects; i++) {
            m_Subjects[i].Initialize(i, labelFiles[i], subjects[0]);
        }
    }

    void ParticleSystem::RunPreprocessing(std::string outputName) {
        double t0 = 0;
        double t1 = 30;
        double dt = 0.01;

        m_ImageContext.ComputeIntersection();

        ParticleSubjectArray subjects(1);
        subjects[0].m_Name = "Intersection";
        subjects[0].NewParticles(m_NumParticlesPerSubject);
        subjects[0].InitializeRandomPoints(m_ImageContext.GetIntersection());

        cout << "Distance map generation ..." << endl;
        ParticleConstraint initialConstraint;
        LabelVector labels;
        labels.push_back(m_ImageContext.GetIntersection());
        initialConstraint.SetImageList(labels);
        cout << "Distance map generation ... done" << endl;

        int k = 0;
        char trackName[128];
        for (double t = t0; t < t1; t += dt) {
            cout << "Preprocessing Time: " << t << endl;
            InternalForce internalForce;
            internalForce.ComputeForce(subjects);
            initialConstraint.ApplyConstraint(subjects);
            UpdateSystem(subjects, dt);
//            sprintf(trackName, "preprocessing_2_%04d.txt", ++k);
//            SaveSystem(trackName, subjects);
        }
        SaveSystem(outputName, subjects);
    }

    void ParticleSystem::Run() {
        const double t0 = 0;
        const double t1 = 30;
        const double dt = 0.01;

        cout << "Distance map generation ..." << endl;
        ParticleConstraint constraint;
        constraint.SetImageList(m_ImageContext.GetLabelVector());
        cout << "Distance map generation ... done" << endl;

        char trackName[128];
        int k = 0;
        for (double t = t0; t < t1; t += dt) {
            cout << "Processing Time: " << t << endl;
            InternalForce internalForce;
            EnsembleForce ensembleForce;
            ensembleForce.SetImageContext(&m_ImageContext);
            internalForce.ComputeForce(m_Subjects);
            ensembleForce.ComputeForce(m_Subjects);
            constraint.ApplyConstraint(m_Subjects);
            UpdateSystem(m_Subjects, dt);
            sprintf(trackName, m_TrackingOutputPattern.c_str(), ++k);
            SaveSystem(trackName, m_Subjects);
        }
    }

    void ParticleSystem::PrepareSystem(ParticleSubjectArray& subjects) {
        if (subjects.size() < 1) {
            return;
        }
        const int nSubjects = subjects.size();
        const int nPoints = subjects[0].m_Particles.size();
        for (int i = 0; i < nSubjects; i++) {
            for (int j = 0; j < nPoints; j++) {
                fordim(k) {
                    subjects[i][j].f[k] = 0;
                }
            }
        }
    }

    void ParticleSystem::UpdateSystem(ParticleSubjectArray& subjects, double dt) {
        const int nSubjects = subjects.size();
        const int nPoints = subjects[0].m_Particles.size();
        for (int i = 0; i < nSubjects; i++) {
            ParticleSubject& iShape = subjects[i];
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

    void ParticleSystem::LoadSystem(std::string filename, int cmd) {
        LoadSystem(filename, m_Subjects, cmd);
    }

    void ParticleSystem::LoadSystem(std::string filename, ParticleSubjectArray& subjects, int cmd) {
        using namespace std;
        int nSubjects = 0;
        int nParticles = 0;

        m_ImageContext.Clear();
        m_TrackingOutputPattern = "/data/Particles/Output/tracking_%04d.txt";

        char cbuf[256];
        ifstream in(filename.c_str());
        if (in.is_open()) {
            string name;
            int value;

            while (in.good())
            {
                in.getline(cbuf, sizeof(cbuf));
                if (in.good()) {
                    stringstream ss(cbuf);
                    ss >> name >> value;
                    if (name == "TrackingOutputPattern") {
                        char buf[128];
                        in.getline(buf, sizeof(buf));
                        if (in.good()) {
                            m_TrackingOutputPattern = buf;
                        }
                    } else if (name == "ParticleDimension:") {
                        if (value != __Dim) {
                            throw "Particle Dimension Mismatch!";
                        }
                    } else if (name == "NumParticlesPerSubject:") {
                        m_NumParticlesPerSubject = value;
                    } else if (name == "LabelImages:") {
                        for (int i = 0; i < value; i++) {
                            char buf[128];
                            in.getline(buf, sizeof(buf));
                            if (in.good()) {
                                if (cmd != 0) {
                                    m_ImageContext.LoadLabel(buf);
                                }
                            }
                            //cout << in.exceptions() << endl;
                        }
                    } else if (name == "IntensityImages:") {
                        for (int i = 0; i < value; i++) {
                            char buf[128];
                            in.getline(buf, sizeof(buf));
                            if (cmd != 0) {
                                m_ImageContext.LoadDoubleImage(buf);
                            }
                        }
                    } else if (name == "Subjects:") {
                        nSubjects = value;
                        subjects.clear();
                        subjects.resize(nSubjects);
                        for (int i = 0; i < nSubjects; i++) {
                            char buf[128];
                            in.getline(buf, sizeof(buf));
                            subjects[i].m_Name = buf;
                        }
                    } else if (name == "Particles:") {
                        nParticles = value;
                        for (int i = 0; i < nSubjects; i++) {
                            subjects[i].NewParticles(nParticles);
                            for (int j = 0; j < nParticles; j++) {
                                Particle& p = subjects[i][j];
                                char buf[128];
                                in.getline(buf, sizeof(buf));
                                if (in.good())
                                {
                                    stringstream ss(buf);
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
            }
        }
    }

    void ParticleSystem::SaveSystem(std::string filename, ParticleSubjectArray& subjects) {
        using namespace std;
        ofstream out(filename.c_str());

        out << "ParticleDimension: " << __Dim << endl;
        out << "NumParticlesPerSubject: " << m_NumParticlesPerSubject << endl;
        StringVector& labels = m_ImageContext.GetFileNames();
        if (labels.size() > 0) {
            out << "LabelImages: " << labels.size() << endl;
            for (int i = 0; i < labels.size(); i++) {
                out << labels[i] << endl;
            }
        }

        StringVector& images = m_ImageContext.GetDoubleImageFileNames();
        if (images.size() > 0) {
            out << "IntensityImages: " << images.size() << endl;
            for (int i = 0; i < images.size(); i++) {
                out << images[i] << endl;
            }
        }

        const int nSubjects = subjects.size();
        const int nParticles = subjects[0].GetNumberOfPoints();
        out << "Subjects: " << subjects.size() << endl;
        for (int i = 0; i < nSubjects; i++) {
            out << subjects[i].m_Name << endl;
        }
        out << "Particles: " << nParticles << endl;
        for (int i = 0; i < nSubjects; i++) {
            for (int j = 0; j < nParticles; j++) {
                Particle& p = subjects[i][j];
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

    void ImageContext::Clear() {
        m_FileNames.clear();
        m_LabelImages.clear();
        m_DistanceMaps.clear();
    }

    void ImageContext::LoadLabel(std::string filename) {
        itkcmds::itkImageIO<LabelImage> io;
        LabelImage::Pointer image = io.ReadImageT(filename.c_str());
        m_LabelImages.push_back(image);

        // set default spacing to 1 to match index and physical coordinate space
        LabelImage::SpacingType defaultSpacing;
        defaultSpacing.Fill(1);
        m_LabelImages.back()->SetSpacing(defaultSpacing);

        m_FileNames.push_back(filename);
    }

    void ImageContext::LoadDoubleImage(std::string filename) {
        itkcmds::itkImageIO<DoubleImage> io;
        DoubleImage::Pointer image = io.ReadImageT(filename.c_str());
        m_Images.push_back(image);

        // set default spacing to 1 to match index and physical coordinate space
        DoubleImage::SpacingType defaultSpacing;
        defaultSpacing.Fill(1);
        m_Images.back()->SetSpacing(defaultSpacing);
        
        m_DoubleImageFileNames.push_back(filename);
    }

    LabelImage::Pointer ImageContext::GetLabel(int j) {
        return m_LabelImages[j];
    }

    DoubleImage::Pointer ImageContext::GetDoubleImage(int j) {
        return m_Images[j];
    }

    LabelImage::Pointer ImageContext::GetIntersection() {
        return m_Intersection;
    }

    void ImageContext::ComputeIntersection() {
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
            const int pixelThreshold = 1;
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
    
    StringVector& ImageContext::GetDoubleImageFileNames() {
        return m_DoubleImageFileNames;
    }
    
    StringVector& ImageContext::GetFileNames() {
        return m_FileNames;
    }
    
    LabelVector& ImageContext::GetLabelVector() {
        return m_LabelImages;
    }
    
    DoubleImageVector& ImageContext::GetDoubleImageVector() {
        return m_Images;
    }
}
