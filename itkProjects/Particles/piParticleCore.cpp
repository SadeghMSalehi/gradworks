#include "iostream"
#include "sstream"
#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/timer.hpp>

#include "piImageIO.h"
#include "piParticleCore.h"
#include "piParticleBSpline.h"
#include "piParticleForces.h"
#include "piParticleCollision.h"
#include "piImageProcessing.h"
#include "vtkMatrix4x4.h"
#include "vtkPolyData.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "piVTK.h"
#include "vtkLandmarkTransform.h"

namespace pi {
    ostream& operator<<(ostream& os, const vtkTransformType& t) {
        vtkMatrix4x4* mat = t->GetMatrix();
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                os << mat->GetElement(i, j) << " ";
            }
        }
        return os;
    }

    istream& operator>>(istream& is, vtkTransformType& t) {
        vtkMatrix4x4* mat = t->GetMatrix();
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                double v;
                is >> v;
                mat->SetElement(i, j, v);
            }
        }
        t->SetMatrix(mat);
        return is;
    }
 
    ParticleSubject::ParticleSubject(): m_SubjId(-1) {
    }

    ParticleSubject::ParticleSubject(int subjid, int npoints) : m_SubjId(subjid) {
        NewParticles(npoints);
    }

    ParticleSubject::~ParticleSubject() {

    }

    void ParticleSubject::Clear() {
        m_Particles.clear();
        m_SubjId = -1;
        m_InverseDeformableTransform = m_DeformableTransform = FieldTransformType::Pointer(NULL);
    }

    void ParticleSubject::Zero() {
        const int nPoints = GetNumberOfPoints();
        for (int i = 0; i < nPoints; i++) {
            m_Particles[i].idx = i;
            m_Particles[i].Zero();
        }
    }

    void ParticleSubject::NewParticles(int n) {
        // FIXME: be careful if something gets wrong
        if (n == m_Particles.size()) {
            Zero();
            return;
        }
        m_Particles.resize(n);
        if (pointscopy.GetPointer() == NULL) {
            pointscopy = vtkPointsType::New();
            pointscopy->SetNumberOfPoints(n);
        }
        Zero();
        if (alignment.GetPointer() == NULL) {
            alignment = vtkTransformType::New();
            inverseAlignment = vtkTransformType::New();
            alignment->Identity();
            inverseAlignment->Identity();
        }
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

        if (indexes.size() <= nPoints) {
            cout << "too small intersection points...[" << indexes.size() << "]" << endl;
            return;
        }

        for (int i = 0; i < nPoints; i++) {
            LabelImage::IndexType idx = indexes[i];
            m_Particles[i].idx = i;
            fordim (k) {
                m_Particles[i].x[k] = idx[k];
            }
        }
    }

    void ParticleSubject::Initialize(int subj, std::string name, int nPoints) {
        m_SubjId = subj;
        if (name != "") {
            m_Name = name;
        }
        NewParticles(nPoints);
    }
    
    void ParticleSubject::Initialize(int subj, std::string name, const ParticleSubject& shape) {
        m_SubjId = subj;
        if (name != "") {
            m_Name = name;
        }
        const int nPoints = shape.GetNumberOfPoints();
        NewParticles(shape.m_Particles.size());
        for (int i = 0; i < nPoints; i++) {
            m_Particles[i] = shape.m_Particles[i];
            m_Particles[i].idx = i;
        }
    }

    void ParticleSubject::Initialize(const ParticleArray& array) {
        const int nPoints = array.size();
        NewParticles(nPoints);
        for (int i = 0; i < nPoints; i++) {
            m_Particles[i] = array[i];
            m_Particles[i].idx = i;
            m_Particles[i].subj = m_SubjId;
        }
    }

    void ParticleSubject::SyncPointsCopy() {
        const int npoints = m_Particles.size();
        for (int i = 0; i < npoints; i++) {
            pointscopy->SetPoint(i, m_Particles[i].x);
        }
    }

    //////////////////////////////////////////////////////////////////////////////////
    //
    // Density computation
    //
    // A density of a particle is defined by the inverse of average neighbor distances
    //
    void ParticleSubject::ComputeDensity() {
        const int npoints = m_Particles.size();
        for (int i = 0; i < npoints; i++) {
            DataReal d = 0;
            for (int j = 0; j < npoints; j++) {
                d += sqrt(m_Particles[i].Dist2(m_Particles[j]));
            }
            m_Particles[i].density = npoints / d;
        }
    }

    void ParticleSubject::ComputeAlignment(ParticleSubject& dst, bool useSimilarity) {
        const bool useAlignment = false;
        if (useAlignment) {
            typedef vtkSmartPointer<vtkLandmarkTransform> vtkLandmarkTransformType;
            dst.SyncPointsCopy();
            SyncPointsCopy();
            vtkLandmarkTransformType landmarkTransform = vtkLandmarkTransformType::New();
            if (useSimilarity) {
                landmarkTransform->SetModeToSimilarity();
            }
            landmarkTransform->SetSourceLandmarks(pointscopy);
            landmarkTransform->SetTargetLandmarks(dst.pointscopy);
            landmarkTransform->Update();

            alignment->SetMatrix(landmarkTransform->GetMatrix());
            inverseAlignment->SetMatrix(alignment->GetMatrix());
            inverseAlignment->Inverse();
        }
    }
    
    void ParticleSubject::AlignmentTransformX2Y() {
        const int npoints = m_Particles.size();
        for (int i = 0; i < npoints; i++) {
            Particle& p = m_Particles[i];
            alignment->TransformPoint(p.x, p.y);
        }
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

    void ParticleSubject::TransformX2X(TransformType* transform) {
        const int nPoints = GetNumberOfPoints();
        if (transform != NULL) {
            TransformType::InputPointType inputPoint;
            for (int i = 0; i < nPoints; i++) {
                fordim (j) {
                    inputPoint[j] = m_Particles[i].x[j];
                }
                TransformType::OutputPointType outputPoint = transform->TransformPoint(inputPoint);
                fordim (j) {
                    m_Particles[i].x[j] = outputPoint[j];
                }
            }
        }
    }

    void ParticleSubject::TransformY2Y(TransformType* transform) {
        const int nPoints = GetNumberOfPoints();
        if (transform != NULL) {
            for (int i = 0; i < nPoints; i++) {
                TransformType::InputPointType inputPoint;
                fordim (j) {
                    inputPoint[j] = m_Particles[i].y[j];
                }
                TransformType::OutputPointType outputPoint = transform->TransformPoint(inputPoint);
                fordim (j) {
                    m_Particles[i].y[j] = outputPoint[j];
                }
            }
        }
    }

    void ParticleSubject::TransformY2Z(TransformType* transform) {
        const int nPoints = GetNumberOfPoints();
        if (transform != NULL) {
            for (int i = 0; i < nPoints; i++) {
                TransformType::InputPointType inputPoint;
                fordim (j) {
                    inputPoint[j] = m_Particles[i].y[j];
                }
                TransformType::OutputPointType outputPoint = transform->TransformPoint(inputPoint);
                fordim (j) {
                    m_Particles[i].z[j] = outputPoint[j];
                }
            }
        } else {
            for (int i = 0; i < nPoints; i++) {
                fordim (k) {
                    m_Particles[i].z[k] = m_Particles[i].y[k];
                }
            }
        }
    }

    void ParticleSubject::TransformZ2Y(TransformType* transform) {
        const int nPoints = GetNumberOfPoints();
        if (transform != NULL) {
            for (int i = 0; i < nPoints; i++) {
                TransformType::InputPointType inputPoint;
                fordim (j) {
                    inputPoint[j] = m_Particles[i].z[j];
                }
                TransformType::OutputPointType outputPoint = transform->TransformPoint(inputPoint);
                fordim (j) {
                    m_Particles[i].y[j] = outputPoint[j];
                }
            }
        } else {
            for (int i = 0; i < nPoints; i++) {
                fordim (k) {
                    m_Particles[i].y[k] = m_Particles[i].z[k];
                }
            }
        }
    }

    void ParticleSubject::ReadParticlePositions(std::istream& is, int nPoints) {
        if (m_Particles.size() != nPoints) {
            m_Particles.resize(nPoints);
        }
        for (int i = 0; i < nPoints; i++) {
            m_Particles[i].idx = i;
            m_Particles[i].t = 0;
            m_Particles[i].subj = m_SubjId;
            for4(k) {
                is >> m_Particles[i].x[k];
            }
        }
    }
    
    void ParticleSubject::WriteParticlePositions(std::ostream& os) {
        for (int i = 0; i < m_Particles.size(); i++) {
            for4(k) {
                os << m_Particles[i].x[k] << "  ";
            }
            os << endl;
        }
    }
    
    void ParticleSubject::ReadParticles(std::istream& is, int nPoints) {
        if (m_Particles.size() != nPoints) {
            m_Particles.resize(nPoints);
        }
        
        for (int i = 0; i < nPoints; i++) {
            m_Particles[i].idx = i;
            m_Particles[i].t = 0;
            m_Particles[i].subj = m_SubjId;
            is >> m_Particles[i];
        }
    }
    
    void ParticleSubject::WriteParticles(std::ostream& os) {
        for (int i = 0; i < m_Particles.size(); i++) {
            os << m_Particles[i] << endl;
        }
    }

    void ParticleSubject::ReadAlignment(std::istream& is) {

    }

    void ParticleSubject::WriteAlignment(std::ostream& os) {

    }

    int ImageContext::Count() {
        return m_Images.size();
    }

    void ImageContext::Clear() {
        m_FileNames.clear();
        m_LabelImages.clear();
        m_DistanceMaps.clear();
    }


    void ImageContext::LoadLabel(std::string filename) {
        ImageIO<LabelImage> io;
        LabelImage::Pointer image = io.ReadCastedImage(filename.c_str());
        if (image.IsNull()) {
            return;
        }
        m_LabelImages.push_back(image);

        // set default spacing to 1 to match index and physical coordinate space
        LabelImage::SpacingType defaultSpacing;
        defaultSpacing.Fill(1);
        m_LabelImages.back()->SetSpacing(defaultSpacing);

        m_FileNames.push_back(filename);
    }

    void ImageContext::LoadRealImage(std::string filename) {
        ImageIO<RealImage> io;
        RealImage::Pointer image = io.ReadCastedImage(filename.c_str());

        ImageProcessing proc;
        image = proc.NormalizeIntensity(image, LabelImage::Pointer(), 0.01);

        // set default spacing to 1 to match index and physical coordinate space
        RealImage::SpacingType defaultSpacing;
        defaultSpacing.Fill(1);
        image->SetSpacing(defaultSpacing);
        
        m_RealImageFileNames.push_back(filename);
        m_Images.push_back(image);
    }

    LabelImage::Pointer ImageContext::GetLabel(int j) {
        return m_LabelImages[j];
    }

    RealImage::Pointer ImageContext::GetRealImage(int j) {
        return m_Images[j];
    }

    LabelImage::Pointer ImageContext::GetIntersection() {
        return m_Intersection;
    }
    
    void ImageContext::SetIntersection(LabelImage::Pointer intersection) {
        m_Intersection = intersection;
    }

    int ImageContext::ComputeIntersection() {
        ImageIO<LabelImage> io;
        LabelImage::Pointer intersection = io.CopyImage(m_LabelImages[0]);
        intersection->FillBuffer(0);

        LabelImage::RegionType region = intersection->GetBufferedRegion();

        // set as member variable to reuse
        m_Intersection = intersection;

        // compute intersection by looping over region
        std::vector<LabelImage::IndexType> indexes;
        LabelImageIteratorType iter(intersection, region);
        int nIntersectionPixels = 0;
        for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
            LabelImage::IndexType idx = iter.GetIndex();
            LabelImage::PixelType pixel = 1;
            const int pixelThreshold = 1;
            for (int i = 0; i < m_LabelImages.size(); i++) {
                if (m_LabelImages[i]->GetPixel(idx) < pixelThreshold) {
                    pixel = 0;
                    break;
                }
            }
            if (pixel > 0) {
                nIntersectionPixels++;
                intersection->SetPixel(idx, 255);
            }
        }
        if (m_IntersectionOutput != "") {
            io.WriteImage(m_IntersectionOutput.c_str(), intersection);
        }
        return nIntersectionPixels;
    }
    
    StringVector& ImageContext::GetRealImageFileNames() {
        return m_RealImageFileNames;
    }
    
    StringVector& ImageContext::GetFileNames() {
        return m_FileNames;
    }
    
    LabelVector& ImageContext::GetLabelVector() {
        return m_LabelImages;
    }
    
    RealImageVector& ImageContext::GetRealImageVector() {
        return m_Images;
    }

    ParticleSystem::ParticleSystem()  {
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

    void ParticleSystem::InitializeSystem(int nsubj, int nparticles) {
        m_Subjects.clear();
        m_Subjects.resize(nsubj);
        for (int i = 0; i < m_Subjects.size(); i++) {
            m_Subjects[i].Initialize(nsubj, "", nparticles);
        }
    }
    
    void ParticleSystem::InitializeSystem(Options& options) {
        m_Options = &options;
        m_Subjects.clear();
        m_Subjects.resize(options.GetStringVector("Subjects:").size());
        for (int i = 0; i < m_Subjects.size(); i++) {
            m_Subjects[i].Initialize(i, options.GetStringVectorValue("Subjects:", i), options.GetInt("NumberOfParticles:", 0));
        }
    }

    // load images for adaptive sampling, or create them
    // adaptive sampling can be turned on by option '+adaptve_sampling'
    void ParticleSystem::LoadKappaImages(Options& options, ImageContext* context) {
        StringVector& kappaNames = options.GetStringVector("KappaImageCache:");
        if (kappaNames.size() != m_Subjects.size()) {
            cout << "Kappa images and subjects are different set (KappaImageCache: missing?)" << endl;
            return;
        }

        DataReal sigma = options.GetReal("AdaptiveSamplingBlurSigma:", 1);
        DataReal maxKappa = options.GetReal("AdaptiveSamplingMaxKappa:", 2);

        ImageIO<RealImage> io;
        for (int i = 0; i < kappaNames.size(); i++) {
            if (io.FileExists(kappaNames[i].c_str())) {
                m_Subjects[i].kappaImage = io.ReadCastedImage(kappaNames[i].c_str());
            } else {
                RealImage::Pointer realImg = context->GetRealImage(i);
                if (realImg.IsNotNull()) {
                    ImageProcessing proc;
                    RealImage::Pointer gradImg = proc.ComputeGaussianGradientMagnitude(realImg, sigma);
                    RealImage::Pointer kappaImg = proc.RescaleIntensity<RealImage>(gradImg, 1.0, maxKappa);
                    m_Subjects[i].kappaImage = kappaImg;
                    io.WriteImage(kappaNames[i].c_str(), kappaImg);
                } else {
                    cout << "Real image is null for subject: " << m_Subjects[i].m_Name << endl;
                }
            }
            if (m_Subjects[i].kappaImage.IsNotNull()) {
                m_Subjects[i].kappaSampler = LinearImageInterpolatorType::New();
                m_Subjects[i].kappaSampler->SetInputImage(m_Subjects[i].kappaImage);
            }
        }
    }

    ParticleSubject& ParticleSystem::GetInitialSubject() {
        return m_InitialSubject;
    }

    ParticleSubject& ParticleSystem::InitializeMean() {
        ComputeXMeanSubject();
        m_MeanSubject.TransformX2Y(NULL);
        m_MeanSubject.TransformY2Z(NULL);
        return m_MeanSubject;
    }

    ParticleSubject& ParticleSystem::ComputeXMeanSubject() {
        const int nPoints = m_Subjects[0].GetNumberOfPoints();
        const int nSubjects = m_Subjects.size();

        m_MeanSubject.m_SubjId = -1;
        m_MeanSubject.NewParticles(nPoints);

        // for every dimension k
        fordim(k) {
            // for every point i
            for (int i = 0; i < nPoints; i++) {
                // sum over all subject j
                for (int j = 0; j < nSubjects; j++) {
                    m_MeanSubject[i].x[k] += m_Subjects[j][i].x[k];
                }
                m_MeanSubject[i].x[k] /= nSubjects;
            }
        }
        return m_MeanSubject;
    }

    ParticleSubject& ParticleSystem::ComputeYMeanSubject() {
        const int nPoints = m_Subjects[0].GetNumberOfPoints();
        const int nSubjects = m_Subjects.size();

        m_MeanSubject.m_SubjId = -1;
        m_MeanSubject.NewParticles(nPoints);

        // for every dimension k
        fordim(k) {
            // for every point i
            for (int i = 0; i < nPoints; i++) {
                // sum over all subject j
                for (int j = 0; j < nSubjects; j++) {
                    m_MeanSubject[i].y[k] += m_Subjects[j][i].y[k];
                }
                m_MeanSubject[i].y[k] /= nSubjects;
            }
        }
        return m_MeanSubject;
    }

    ParticleSubject& ParticleSystem::ComputeZMeanSubject() {
        const int nPoints = m_Subjects[0].GetNumberOfPoints();
        const int nSubjects = m_Subjects.size();

        m_MeanSubject.m_SubjId = -1;
        m_MeanSubject.NewParticles(nPoints);

        // for every dimension k
        fordim(k) {
            // for every point i
            for (int i = 0; i < nPoints; i++) {
                // sum over all subject j
                for (int j = 0; j < nSubjects; j++) {
                    m_MeanSubject[i].z[k] += m_Subjects[j][i].z[k];
                }
                m_MeanSubject[i].z[k] /= nSubjects;
            }
        }
        return m_MeanSubject;
    }

    ParticleSubject& ParticleSystem::GetMeanSubject() {
        return m_MeanSubject;
    }

    ParticleSubjectArray& ParticleSystem::GetSubjects() {
        return m_Subjects;
    }

    void export2vtk(ParticleSubject& sub, const char* vtkname, int field) {
        int np = sub.m_Particles.size();
        vtkPolyData* vtk = vtkPolyData::New();
        vtkPoints* points = vtkPoints::New();
        vtkFloatArray* f = vtkFloatArray::New();
        f->SetName("F");
        f->SetNumberOfComponents(3);
        f->SetNumberOfTuples(np);
        vtkFloatArray* v = vtkFloatArray::New();
        v->SetName("V");
        v->SetNumberOfComponents(3);
        v->SetNumberOfTuples(np);
        vtkFloatArray* w = vtkFloatArray::New();
        w->SetName("W");
        w->SetNumberOfComponents(3);
        w->SetNumberOfTuples(np);
        vtk->SetPoints(points);

        bool useV = (field&2) > 0;
        if (useV) {
            vtk->GetPointData()->AddArray(v);
        }
        bool useF = (field&4) > 0;
        if (useF) {
            vtk->GetPointData()->AddArray(f);
        }
        const bool useW = (field&8) > 0;
        if (useW) {
            vtk->GetPointData()->AddArray(w);
        }
        points->SetNumberOfPoints(np);
        // field \in [0,15]
        for (int j = 0; j < np; j++) {
            Particle& p = sub[j];
            if ((field & 1) == 0) {
                points->SetPoint(j, p.x);
            }
            if ((field & 1) > 0) {
                points->SetPoint(j, p.y);
            }
            if (useV) {
                v->SetTuple(j, p.v);
            }
            if (useF) {
                f->SetTuple(j, p.f);
            }
            if (useW) {
                w->SetTuple(j, p.w);
            }
        }
        pivtk::vtk_write_polydata(vtkname, vtk);
        v->Delete();
        f->Delete();
        w->Delete();
        points->Delete();
        vtk->Delete();
    }


    vtkPolyData* convert2vtk(ParticleArray& parray) {
        vtkPolyData* out = vtkPolyData::New();
        vtkPoints* points = vtkPoints::New();
        points->SetNumberOfPoints(parray.size());
        for (int i = 0; i < parray.size(); i++) {
            points->SetPoint(i, parray[i].x);
        }
        out->SetPoints(points);
        return out;
    }
}
