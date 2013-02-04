#include "myParticleCore.h"

#include "myParticleConstraint.h"
#include "myParticleBSpline.h"
#include "iostream"
#include "sstream"

#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/timer.hpp>

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

    void Particle::AddForce(const double* ff, double alpha) {
        fordim(i) {
            f[i] += (alpha * ff[i]);
        }
    }

    void Particle::SubForce(const double* ff, double alpha) {
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

    void ParticleSubject::Clear() {
        m_Particles.clear();
        m_SubjId = -1;
        m_InverseTransform = m_Transform = FieldTransformType::Pointer(NULL);
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
        if (name != "") {
            m_Name = name;
        }
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

    void ParticleSlice::Update(ParticleSubjectArray& subjects, LabelImage::Pointer labelImage) {
        LabelImage::SizeType sz = labelImage->GetBufferedRegion().GetSize();

        const int nSubj = subjects.size();
        const int nPoints = subjects[0].GetNumberOfPoints();

        m_ParticlePointerMatrix.resize(sz[m_SliceDim], nSubj);
        for (int i = 0; i < nSubj; i++) {
            for (int j = 0; j < nPoints; j++) {
                for (int k = 0; k < sz[m_SliceDim]; k++) {
                    if (subjects[i][j].x[m_SliceDim] >= k && subjects[i][j].x[m_SliceDim] < k) {
                        m_ParticlePointerMatrix(k, i).push_back(&subjects[i][j]);
                        break;
                    }
                }
            }
        }
    }

    const ParticleSlice::ParticlePointerVector&  ParticleSlice::Get(int slice, int subj) {
        return m_ParticlePointerMatrix(slice, subj);
    }

    ParticleSystem::ParticleSystem() : m_NumParticlesPerSubject(300) {
        m_times[0] = 0;
        m_times[1] = 0.01;
        m_times[2] = 10;

        m_timesPreprocessing[0] = 0;
        m_timesPreprocessing[1] = 0.01;
        m_timesPreprocessing[2] = 10;
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

    /*
    void ParticleSystem::LoadPreprocessing(std::string filename) {
        ParticleSubjectArray subjects;
        if (filename != "") {
            LoadSystem(filename, m_Initial, 1);
            StringVector& labelFiles = m_ImageContext.GetFileNames();
            const int nSubjects = labelFiles.size();
            m_Subjects.resize(nSubjects);
            for (int i = 0; i < nSubjects; i++) {
                m_Subjects[i].Initialize(i, labelFiles[i], subjects[0]);
            }
        }
    }
    */

    void ParticleSystem::RunPreprocessing() {
        if (m_Initial.size() == 0) {
            double t0 = m_timesPreprocessing[0];
            double dt = m_timesPreprocessing[1];
            double t1 = m_timesPreprocessing[2];

            m_ImageContext.ComputeIntersection();

            m_Initial.resize(1);
            m_Initial[0].m_Name = "Intersection";
            m_Initial[0].NewParticles(m_NumParticlesPerSubject);
            m_Initial[0].InitializeRandomPoints(m_ImageContext.GetIntersection());

            cout << "Distance map generation ..." << endl;
            ParticleConstraint initialConstraint;
            LabelVector labels;
            labels.push_back(m_ImageContext.GetIntersection());
            initialConstraint.SetImageList(labels);
            cout << "Distance map generation ... done" << endl;

            int k = 0;
            char trackName[128];
            for (double t = t0; t < t1; t += dt) {
                boost::timer timer;
                InternalForce internalForce;
                internalForce.ComputeForce(m_Initial);
                initialConstraint.ApplyConstraint(m_Initial);
                UpdateSystem(m_Initial, dt);
                //            sprintf(trackName, "preprocessing_2_%04d.txt", ++k);
                //            SaveSystem(trackName, subjects);
                cout << "Preprocessing Time: " << t << "; Elapsed Time: " << timer.elapsed() << " secs" << endl;
            }
        }
    }

    void ParticleSystem::Run() {
        const double t0 = m_times[0];
        const double dt = m_times[1];
        const double t1 = m_times[2];


        cout << "Distance map generation ..." << endl;
        ParticleConstraint constraint;
        constraint.SetImageList(m_ImageContext.GetLabelVector());
        cout << "Distance map generation ... done" << endl;

        char trackName[128];
        int k = 0;
        boost::timer timer;
        for (double t = t0; t < t1; t += dt) {
            timer.restart();
            InternalForce internalForce;
            EnsembleForce ensembleForce;

            internalForce.ComputeForce(m_Subjects);
            ensembleForce.SetImageContext(&m_ImageContext);
            ensembleForce.ComputeEnsembleForce(m_Subjects);
            ensembleForce.ComputeImageForce(m_Subjects);

            constraint.ApplyConstraint(m_Subjects);
            UpdateSystem(m_Subjects, dt);
            if (m_TrackingOutputPattern != "") {
                sprintf(trackName, m_TrackingOutputPattern.c_str(), ++k);
                SaveSystem(trackName);
            }
            cout << "Processing Time: " << t << "; Elapsed " << timer.elapsed() << " secs" << endl;
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

    bool ParticleSystem::LoadSystem(std::string filename) {
        using namespace std;
        int nSubjects = 0;

        m_ImageContext.Clear();
        m_TrackingOutputPattern = "";

        m_Subjects.clear();
        m_Initial.clear();

        char cbuf[256];
        ifstream in(filename.c_str());
        if (in.is_open()) {
            string name;
            while (in.good()) {
                in.getline(cbuf, sizeof(cbuf));
                cout << cbuf << endl;
                if (in.good()) {
                    stringstream ss(cbuf);
                    ss >> name;
                    if (name[0] == '#') {
                        continue;
                    } else if (name == "TrackingOutputPattern:") {
                        ss >> m_TrackingOutputPattern;
                    } else if (name == "ParticleDimension:") {
                        int value;
                        ss >> value;
                        if (value != __Dim) {
                            cout << "Particle Domain Mismatch!" << endl;
                            throw "Particle Dimension Mismatch!";
                        }
                    } else if (name == "NumParticlesPerSubject:") {
                        ss >> m_NumParticlesPerSubject;
                    } else if (name == "TimeRange:") {
                        ss >> m_times[0] >> m_times[1] >> m_times[2];
                    } else if (name == "PreprocessingTimeRange:") {
                        ss >> m_timesPreprocessing[0] >> m_timesPreprocessing[1] >> m_timesPreprocessing[2];
                    } else if (name == "LabelImages:") {
                        int value;
                        ss >> value;
                        for (int i = 0; i < value; i++) {
                            char buf[128];
                            in.getline(buf, sizeof(buf));
                            if (in.good()) {
                                try {
                                    m_ImageContext.LoadLabel(buf);
                                } catch (itk::ExceptionObject& ex) {
                                    ex.Print(cout);
                                    exit(0);
                                }
                            }
                            //cout << in.exceptions() << endl;
                        }
                    } else if (name == "IntensityImages:") {
                        int value;
                        ss >> value;
                        for (int i = 0; i < value; i++) {
                            char buf[128];
                            in.getline(buf, sizeof(buf));
                            m_ImageContext.LoadDoubleImage(buf);
                        }
                    } else if (name == "Subjects:") {
                        int value;
                        ss >> value;
                        nSubjects = value;

                        m_Subjects.clear();
                        m_Subjects.resize(nSubjects);
                        for (int i = 0; i < nSubjects; i++) {
                            char buf[128];
                            in.getline(buf, sizeof(buf));
                            m_Subjects[i].m_Name = buf;
                        }
                    } else if (name == "InitialParticles:") {
                        m_Initial.resize(1);
                        for (int i = 0; i < m_Initial.size(); i++) {
                            m_Initial[i].NewParticles(m_NumParticlesPerSubject);
                            for (int j = 0; j < m_NumParticlesPerSubject; j++) {
                                Particle& p = m_Initial[i][j];
                                char buf[128];
                                in.getline(buf, sizeof(buf));
                                if (in.good())
                                {
                                    stringstream ss(buf);
                                    for4(k) { ss >> p.x[k]; }
//                                    for4(k) { ss >> p.y[k]; }
//                                    for4(k) { ss >> p.v[k]; }
//                                    for4(k) { ss >> p.f[k]; }
//                                    ss >> p.density;
//                                    ss >> p.pressure;
                                }
                            }
                        }

                    } else if (name == "Particles:") {
                        int subjId;
                        ss >> subjId;
                        int nParticles;
                        ss >> nParticles;
                        if (nParticles > 0) {
                            cout << "Reading " << nParticles << " particles ..." << endl;
                            if (m_Subjects[subjId].GetNumberOfPoints() != nParticles) {
                                m_Subjects[subjId].NewParticles(m_NumParticlesPerSubject);
                            }
                            for (int j = 0; j < nParticles; j++) {
                                Particle& p = m_Subjects[subjId][j];
                                char buf[128];
                                in.getline(buf, sizeof(buf));
                                if (in.good()) {
                                    stringstream ss(buf);
                                    for4(k) { ss >> p.x[k]; }
                                    for4(k) { ss >> p.y[k]; }
                                    for4(k) { ss >> p.v[k]; }
                                    for4(k) { ss >> p.f[k]; }
                                    ss >> p.density;
                                    ss >> p.pressure;
                                }
                            }
                        } else if (nParticles == 0) {
                            cout << "Using " << m_Initial[0].GetNumberOfPoints() << " initial particles ..." << endl;
                            m_Subjects[subjId].Initialize(subjId, "", m_Initial[0]);
                        }
                    }
                }
            }
        } else {
            return false;
        }
        return true;
    }

    void ParticleSystem::SaveSystem(std::string filename) {
        using namespace std;
        ofstream out(filename.c_str());

        out << "ParticleDimension: " << __Dim << endl;
        out << "NumParticlesPerSubject: " << m_NumParticlesPerSubject << endl;

        // Time Range
        out << "PreprocessingTimeRange: " << m_timesPreprocessing[0] << " " << m_timesPreprocessing[1] << " " << m_timesPreprocessing[2] << endl;

        out << "TimeRange: " << m_times[0] << " " << m_times[1] << " " << m_times[2] << endl;

        // Tracking output pattern
        if (m_TrackingOutputPattern != "") {
            out << "TrackingOutputPattern: " << m_TrackingOutputPattern << endl;
        }

        // Intersection file name
        if (m_ImageContext.m_IntersectionOutput != "") {
            out << "IntersectionImage: " << m_ImageContext.m_IntersectionOutput << endl;
        }

        const int nSubjects = m_Subjects.size();
        out << "Subjects: " << nSubjects << endl;
        for (int i = 0; i < nSubjects; i++) {
            out << m_Subjects[i].m_Name << endl;
        }
        
        // Label Images
        StringVector& labels = m_ImageContext.GetFileNames();
        if (labels.size() > 0) {
            out << "LabelImages: " << labels.size() << endl;
            for (int i = 0; i < labels.size(); i++) {
                out << labels[i] << endl;
            }
        }

        // Intensity Images
        StringVector& images = m_ImageContext.GetDoubleImageFileNames();
        if (images.size() > 0) {
            out << "IntensityImages: " << images.size() << endl;
            for (int i = 0; i < images.size(); i++) {
                out << images[i] << endl;
            }
        }

        // Initial Particles
        out << "InitialParticles: " << m_NumParticlesPerSubject << endl;
        for (int i = 0; i < m_Initial.size(); i++) {
            for (int j = 0; j < m_NumParticlesPerSubject; j++) {
                Particle& p = m_Initial[i][j];
                for4(k) { out << p.x[k] << " "; }
//                for4(k) { out << p.y[k] << " "; }
//                for4(k) { out << p.v[k] << " "; }
//                for4(k) { out << p.f[k] << " "; }
//                out << p.density << " ";
//                out << p.pressure << " ";
                out << endl;
            }
        }

        // Particles
//        cout << "Subjects: " << m_Subjects.size() << endl;
        for (int i = 0; i < m_Subjects.size(); i++) {
            int nParticles = m_Subjects[i].GetNumberOfPoints();
            out << "Particles: " << i << " " << nParticles << endl;
            for (int j = 0; j < nParticles; j++) {
                Particle& p = m_Subjects[i][j];
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
        if (m_IntersectionOutput != "") {
            io.WriteImageT(m_IntersectionOutput.c_str(), intersection);
        }
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
