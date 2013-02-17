//
//  myParticleForces.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 1/26/13.
//
//

#include "piParticleCore.h"
#include "piParticleBSpline.h"
#include "piParticleForces.h"
#include "piParticleSystem.h"
#include "itkGradientImageFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkConstNeighborhoodIterator.h"
#include "armadillo"


#define __show2(t,x,m,n) for (int _=0;_<m;_++) { cout << t << " #" << _ << ":"; for (int __=0;__<n;__++) cout << " " << x[__+_*n]; cout << endl; }
#define __show1(x,m) cout << #x << ":"; for (int _=0;_<m;_++) cout << " " << x[_]; cout << endl;
#define __showmatrix(x,m,n) for (int _=0;_<m;_++) { for (int __=0;__<n;__++) cout << x(_,__); cout << endl; }; cout << endl
#define __showcmd(x,m,n,cmd) for (int _=0;_<m;_++) { for (int __=0;__<n;__++) cmd; cout << endl; }; cout << endl

namespace pi {
    typedef itk::GradientImageFilter<DoubleImage> GradientFilterType;
    typedef GradientFilterType::OutputImageType GradientImage;
    typedef GradientFilterType::OutputPixelType GradientPixel;
    typedef itk::ConstNeighborhoodIterator<GradientImage> VectorImageNeighborhoodIteratorType;

    typedef itk::GradientRecursiveGaussianImageFilter<DoubleImage, GradientImage> GaussianGradientFilterType;
    typedef itk::VectorLinearInterpolateImageFunction<GradientPixel> GradientInterpolatorType;
    typedef itk::ConstNeighborhoodIterator<DoubleImage> DoubleImageNeighborhoodIteratorType;
    typedef itk::ConstNeighborhoodIterator<GradientImage> GradientImageNeighborhoodIteratorType;

//    static void ExtractAttributes(DoubleImage::Pointer image, VectorImage::Pointer grad, Particle& par);

    void InternalForce::ComputeForce(ParticleSubject& subj) {
        const int nPoints = subj.m_Particles.size();
        ParticleArray& particles = subj.m_Particles;
#pragma omp parallel for
        for (int i = 0; i < nPoints; i++) {
            Particle& pi = particles[i];
            for (int j = i+1; j < nPoints; j++) {
                Particle& pj = particles[j];
                ComputeForce(pi, pj);
            }
        }
    }
    
    void InternalForce::ComputeForce(ParticleSubjectArray& shapes) {
        const int nSubjects = shapes.size();
        for (int k = 0; k < nSubjects; k++) {
            ComputeForce(shapes[k]);
        }
    }
    
    void InternalForce::ComputeForce(Particle &pi, Particle &pj) {
        const double sigma = 15 * 5;
        const double coeff = M_PI_2 / sigma;
        const bool useSimpleForce = false;
        
        double fi[__Dim] = { 0 }, fj[__Dim] = { 0 };
        double dx[__Dim] = { 0 };

        double rij2 = 0;
        pi.Sub(pj, dx);
        fordim(k) {
            rij2 += (dx[k]*dx[k]);
        }
        const double rij = std::sqrt(rij2);

        if (rij <= sigma) {
            if (useSimpleForce) {
                fordim (k) {
                    fj[k] = fi[k] = -dx[k];
                }
            } else {
                fordim(k) {
                    dx[k] /= rij;
                }
                const double crij = rij * coeff;
                const double sin1crij = std::sin(crij);
                const double sin2crij = sin1crij * sin1crij;
                fordim(k) {
                    fj[k] = fi[k] = (dx[k] * (coeff * (1 - (1 / sin2crij))));
                }
            }
            pi.SubForce(fi);
            pj.AddForce(fj);
        }
    }

    
    void EntropyInternalForce::ComputeForce(ParticleSubject& subj) {
        const int nPoints = subj.GetNumberOfPoints();
        for (int i = 0; i < nPoints; i++) {
            Particle& pi = subj.m_Particles[i];
            // iteration over particles
            // may reduce use symmetric properties
            VNLVector weights(nPoints, 0);
            for (int j = 0; j < nPoints; j++) {
                Particle& pj = subj.m_Particles[j];
                if (i == j) {
                    // there's no self interaction
                    weights[j] = 0;
                } else {
                    const double sigma = 3;
                    double kappa = 1;
                    double sigma2 = sigma * sigma;
                    double cutoff = 15;
                    /*
                     // kappa should use jPos
                     VNLVectorRef jPos(2, &gPos[n][nDim*j]);
                     ContinuousIndexType jIdx;
                     jIdx[0] = jPos[0];
                     jIdx[1] = jPos[1];
                     double kappa = 1;
                     if (useAdaptiveSampling) {
                     kappaIntp->EvaluateAtContinuousIndex(jIdx);
                     kappa *= kappa;
                     }
                     */
                    double dij = sqrt(pi.Dist2(pj));
                    if (dij > cutoff) {
                        weights[j] = 0;
                    } else {
                        weights[j] = exp(-dij*dij*kappa/(sigma2));
                    }
                }
            }
            double sumForce = weights.sum();
            if (sumForce > 0) {
                weights /= sumForce;
            }
            
            // update force for neighboring particles
            for (int j = 0; j < nPoints; j++) {
                if (i == j || weights[j] == 0) {
                    continue;
                }
                Particle& pj = subj.m_Particles[j];
                double weight = weights[j];
                VNLVector xixj(__Dim, 0);
                fordim (k) {
                    xixj[k] = pi.x[k] - pj.x[k];
                }
                xixj.normalize();
                fordim (k) {
                    pi.f[k] += (weight * (pi.x[k] - pj.x[k]));
                }
            }
        }
    }
    
    void EntropyInternalForce::ComputeForce(ParticleSubjectArray& subjs) {
        const int nSubjs = subjs.size();
        for (int n = 0; n < nSubjs; n++) {
            ComputeForce(subjs[n]);
        }
    }

    void EntropyInternalForce::ComputeForce(Particle& a, Particle& b) {

    }


    EnsembleForce::EnsembleForce(double coeff) : m_Coeff(coeff) {

    }

    EnsembleForce::~EnsembleForce() {

    }

    void EnsembleForce::SetImageContext(ImageContext* context) {
        m_ImageContext = context;
    }

    void EnsembleForce::ComputeMeanShape(ParticleSubjectArray& shapes) {
        if (shapes.size() < 1) {
            return;
        }
        
        const int nSubjects = shapes.size();
        const int nPoints = shapes[0].GetNumberOfPoints();

        m_MeanShape.m_SubjId = -1;
        m_MeanShape.NewParticles(shapes[0].GetNumberOfPoints());

        // for every dimension k
        fordim(k) {
            // for every point i
            for (int i = 0; i < nPoints; i++) {
                // sum over all subject j
                for (int j = 0; j < nSubjects; j++) {
                    m_MeanShape[i].x[k] += shapes[j][i].x[k];
                }
                m_MeanShape[i].x[k] /= nSubjects;
            }
        }
    }
    
    void EnsembleForce::ComputeEnsembleForce(ParticleSystem& system) {
        if (system.GetNumberOfSubjects() < 2) {
            return;
        }
        const int nPoints = system.GetNumberOfParticles();
        const int nSubjects = system.GetNumberOfSubjects();

        // Compute offset for origin-centered shapes
        VNLMatrix centers(nSubjects, 4);
        centers.fill(0);
        for (int i = 0; i < nSubjects; i++) {
            ParticleSubject& subject = system[i];
            fordim(k) {
                for (int j = 0; j < nPoints; j++) {
                    Particle& par = subject[j];
                    centers[i][k] += par.x[k];
                }
                centers[i][k] /= nPoints;
            }
            AffineTransformType::Pointer affineTransform = AffineTransformType::New();
            AffineTransformType::OutputVectorType offset;
            fordim(k) {
                offset[k] -= centers[i][k];
            }
            affineTransform->Translate(offset);
            subject.m_AffineTransform = affineTransform;
//            subject.TransformX2X(affineTransform.GetPointer());
        }
        
        // Compute xMeanShape
        ComputeMeanShape(system.GetSubjects());
        for (int i = 0; i < nSubjects; i++) {
            ParticleBSpline bspline;
            bspline.SetReferenceImage(m_ImageContext->GetLabel(i));
            bspline.EstimateTransform(system[i], m_MeanShape);
            FieldTransformType::Pointer deformableTransform = bspline.GetTransform();
            system[i].TransformX2Y(deformableTransform.GetPointer());
            system[i].m_DeformableTransform = deformableTransform;
        }

        // Compute yMeanShape
        fordim(k) {
            for (int j = 0; j < nPoints; j++) {
                for (int i = 0; i < nSubjects; i++) {
                    m_MeanShape[j].y[k] += system[i][j].y[k];
                }
                m_MeanShape[j].y[k] /= nSubjects;
            }
        }

        for (int i = 0; i < nSubjects; i++) {
            ParticleSubject& iSubj = system[i];
            #pragma omp parallel for
            for (int j = 0; j < nPoints; j++) {
                FieldTransformType::InputPointType xPoint;
                FieldTransformType::JacobianType xJac;
                xJac.set_size(__Dim,__Dim);

                double f[__Dim] = { 0, };
                double *x = iSubj.m_Particles[j].x;
                double *y = iSubj.m_Particles[j].y;
                double *my = m_MeanShape[j].y;
                fordim(k) {
                    xPoint[k] = x[k];
                }
                iSubj.m_DeformableTransform->ComputeInverseJacobianWithRespectToPosition(xPoint, xJac);
//                shapes[i].m_Transform->ComputeJacobianWithRespectToPosition(xPoint, xJac);
                fordim(k) {
                    if (__Dim == 3) {
                        f[k] = xJac[0][k]*(y[0]-my[0]) + xJac[1][k]*(y[1]-my[1]) + xJac[2][k]*(y[2]-my[2]);
                    } else if (__Dim == 2) {
                        f[k] = xJac[0][k]*(y[0]-my[0]) + xJac[1][k]*(y[1]-my[1]);
                    }
                }
                iSubj.m_Particles[j].SubForce(f, m_Coeff);
            }
        }

        // recover coordinates of particles
        for (int i = 0; i < nSubjects; i++) {
            ParticleSubject& subject = system[i];
            AffineTransformType::Pointer inverse = AffineTransformType::New();
            subject.m_AffineTransform->GetInverse(inverse.GetPointer());
//            subject.TransformX2X(inverse);
        }
    }



    //

#ifdef DIMENSION3
    template <int N>
    class ParticleAttribute {
    public:
        const static int NATTRS = N*N*N;
        double f[__Dim];
        double F[__Dim];
        double x[NATTRS];
        double y[NATTRS];
        GradientPixel* gptr[NATTRS];
    };
#else
    template <int N>
    class ParticleAttribute {
    public:
        const static int NATTRS = N*N;
        double f[__Dim];
        double F[__Dim];
        double x[NATTRS];
        double y[NATTRS];
        GradientPixel* gptr[NATTRS];
    };
#endif
    typedef ParticleAttribute<3> Attr;
    template <int N>
    ostream& operator<<(ostream& os, ParticleAttribute<N>& attr) {
        int p = os.precision();
        os.precision(2);
        for (int i = 0; i < attr.NATTRS; i++) {
            os << " " << attr.x[i];
        }
        os << ";";
        for (int i = 0; i < attr.NATTRS; i++) {
            os << " " << attr.y[i];
        }
        os << ";";
        for (int i = 0; i < attr.NATTRS; i++) {
            os << " " << attr.gptr[i]->GetElement(0) << "," << attr.gptr[i]->GetElement(1);
        }
        os.precision(p);
        return os;
    }

    IntensityForce::IntensityForce(double coeff) : m_Coeff(coeff) {

    }

    IntensityForce::~IntensityForce() {

    }

    void IntensityForce::SetImageContext(ImageContext* context) {
        m_ImageContext = context;
    }

    void IntensityForce::ComputeIntensityForce(ParticleSystem* system) {
        ParticleSubjectArray& shapes = system->GetSubjects();
        const ParticleSubject& meanSubject = system->GetMeanSubject();
        
        const int nSubj = shapes.size();
        const int nPoints = shapes[0].GetNumberOfPoints();
        const int nRadius = 1;

        DoubleImageVector warpedImages(nSubj);
//        VectorImageVector gradImages(nSubj);
        std::vector<GradientImage::Pointer> gradImages(nSubj);

        typedef boost::numeric::ublas::matrix<Attr> AttrMatrix;
        AttrMatrix attrs(nSubj, nPoints);
        VNLMatrix attrsMean(nPoints, Attr::NATTRS);

        itkcmds::itkImageIO<DoubleImage> io;        
        // first create a warped image into the mean transform space
        // second create a gradient vector image per subject
        // third extract attributes (features)
        for (int i = 0; i < nSubj; i++) {
            ParticleSubject& subject = shapes[i];
            LabelImage::Pointer refImage = m_ImageContext->GetLabel(i);
            
            // warp the real image
            ParticleBSpline bspline;
            bspline.SetReferenceImage(refImage);
            bspline.EstimateTransform(meanSubject, subject);
            FieldTransformType::Pointer deformableTransform = bspline.GetTransform();
            subject.m_InverseDeformableTransform = deformableTransform;
            warpedImages[i] = bspline.WarpImage(m_ImageContext->GetDoubleImage(i));


            if (subject.m_DeformableTransform.IsNull()) {
                ParticleBSpline forwardBspline;
                forwardBspline.SetReferenceImage(refImage);
                forwardBspline.EstimateTransform(subject, meanSubject);
                FieldTransformType::Pointer forwardTransform = bspline.GetTransform();
                subject.TransformX2Y(forwardTransform);
            }

            char warpedname[128];
            sprintf(warpedname, "warped%d_%03d.nrrd", i, system->currentIteration);
            io.WriteImageT(warpedname, warpedImages[i]);

            const bool useGaussianGradient = false;
            if (useGaussianGradient) {
                GaussianGradientFilterType::Pointer grad = GaussianGradientFilterType::New();
                grad->SetInput(warpedImages[i]);
                grad->Update();
                gradImages[i] = grad->GetOutput();
            } else {
                GradientFilterType::Pointer grad = GradientFilterType::New();
                grad->SetInput(warpedImages[i]);
                grad->Update();
                gradImages[i] = grad->GetOutput();
            }

            DoubleImage& warpedImage = *(warpedImages[i].GetPointer());
            GradientImage& gradImage = *(gradImages[i].GetPointer());

            // extract attributes
            DoubleImage::SizeType radius;
            radius.Fill(nRadius);
            DoubleImageNeighborhoodIteratorType iiter(radius, warpedImages[i], warpedImages[i]->GetBufferedRegion());
            DoubleImage::IndexType idx;
            #pragma omp parallel for
            for (int j = 0; j < nPoints; j++) {
                Particle& par = subject.m_Particles[j];
                fordim(k) {
                    idx[k] = round(par.y[k]);
                }
                iiter.SetLocation(idx);
                Attr& jAttr = attrs(i, j);
                for (int k = 0; k < Attr::NATTRS; k++) {
                    IntIndex idx = iiter.GetIndex(k);
                    jAttr.x[k] = warpedImage[idx];
                    jAttr.gptr[k] = &gradImage[idx];
                }
            }
        }
//        
//        for (int j = 0; j < nSubj; j++) {
//            cout << "Attributes #" << j << ": ";
//            for (int k = 0; k < nAttrs; k++) {
//                cout << attrs[k + j * nAttrs] << " ";
//            }
//            cout << endl;
//        }
//        

        // column sum
        #pragma omp parallel for
        for (int i = 0; i < nPoints; i++) {
            for (int j = 0; j < Attr::NATTRS; j++) {
                attrsMean[i][j] = 0;
                for (int s = 0; s < nSubj; s++) {
                    double attr = attrs(s, i).x[j];
                    attrsMean[i][j] += attr;
                }
            }
        }
        attrsMean /= nSubj;
        cout << attrsMean << endl;

        // compute mean differences
        #pragma omp parallel for
        for (int i = 0; i < nSubj; i++) {
            for (int j = 0; j < nPoints; j++) {
                for (int k = 0; k < Attr::NATTRS; k++) {
                    attrs(i,j).y[k] = attrs(i,j).x[k] - attrsMean[j][k];
                }
            }
        }
        __showmatrix(attrs, nSubj, nPoints);

        // compute force direction
        for (int i = 0; i < nSubj; i++) {
            #pragma omp parallel for
            for (int j = 0; j < nPoints; j++) {
                Attr& attr = attrs(i,j);
                fordim (k) {
                    attr.f[k] = 0;
                    for (int l = 0; l < Attr::NATTRS; l++) {
                        GradientPixel& g = *(attr.gptr[l]);
                        double a = attr.y[l];
                        double b = g[k];
                        if (b > 1e-6 || b < -1e-6) {
                            attr.f[k] += (a / g[k]);
                        }
                    }
                }
            }
        }
        __showcmd(attrs, nSubj, nPoints, cout << " " << attrs(_,__).f[0] << "," << attrs(_,__).f[1]);

        // covariance matrix
        

        // compute force at subject space
        for (int i = 0; i < nSubj; i++) {
            FieldTransformType::Pointer fieldTransform = shapes[i].m_InverseDeformableTransform;
            for (int j = 0; j < nPoints; j++) {
                double* forcePtr = attrs(i,j).f;
                double* forceOutPtr = attrs(i,j).F;
                FieldTransformType::InputPointType x;
                fordim(k) {
                    x[k] = shapes[i][j].x[k];
                } 
                FieldTransformType::JacobianType jac;
                jac.set_size(__Dim, __Dim);
                fieldTransform->ComputeInverseJacobianWithRespectToPosition(x, jac);
                fordim(k) {
                    double ff = 0;
                    if (__Dim == 3) {
                        ff = jac[0][k]*forcePtr[k] + jac[1][k]*forcePtr[k] + jac[2][k]*forcePtr[k];
                    } else if (__Dim == 2) {
                        ff = jac[0][k]*forcePtr[k] + jac[1][k]*forcePtr[k];
                    }
                    forceOutPtr[k] = ff;
                }
            }
        }
        __showcmd(attrs, nSubj, nPoints, cout << " " << attrs(_,__).F[0] << "," << attrs(_,__).F[1]);


        for (int i = 0; i < nSubj; i++) {
            ParticleSubject& subj = shapes[i];
            for (int j = 0; j < nPoints; j++) {
                Particle& par = subj[j];
                VNLVector ff(__Dim);
                fordim (k) {
                    ff[k] = attrs(i,j).F[k];
                }
                ff.normalize();
                par.AddForce(ff.data_block(), m_Coeff);
            }
        }
    }
}
