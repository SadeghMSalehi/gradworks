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
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkConstNeighborhoodIterator.h"
#include "armadillo"
#include "boost/numeric/ublas/matrix.hpp"
#include "boost/numeric/ublas/io.hpp"

#define __show2(t,x,m,n) for (int _=0;_<m;_++) { cout << t << " #" << _ << ":"; for (int __=0;__<n;__++) cout << " " << x[__+_*n]; cout << endl; }
#define __show1(x,m) cout << #x << ":"; for (int _=0;_<m;_++) cout << " " << x[_]; cout << endl;
#define __showmatrix(x,m,n) for (int _=0;_<m;_++) { for (int __=0;__<n;__++) cout << x(_,__); cout << endl; }; cout << endl
#define __showcmd(x,m,n,cmd) for (int _=0;_<m;_++) { for (int __=0;__<n;__++) cmd; cout << endl; }; cout << endl

namespace pi {

    typedef itk::ConstNeighborhoodIterator<GradientImage> VectorImageNeighborhoodIteratorType;
    typedef itk::GradientRecursiveGaussianImageFilter<DoubleImage, GradientImage> GaussianGradientFilterType;
    typedef itk::VectorLinearInterpolateImageFunction<GradientPixel> GradientInterpolatorType;
    typedef itk::ConstNeighborhoodIterator<DoubleImage> DoubleImageNeighborhoodIteratorType;
    typedef itk::ConstNeighborhoodIterator<GradientImage> GradientImageNeighborhoodIteratorType;

    ostream& operator<<(ostream& os, Attr& attr) {
        int p = os.precision();
        os.precision(8);
        for (int i = 0; i < attr.NATTRS; i++) {
            os << " " << attr.x[i];
        }
        os << ";";
        for (int i = 0; i < attr.NATTRS; i++) {
            os << " " << attr.y[i];
        }
        os << ";";
        for (int i = 0; i < attr.NATTRS; i++) {
            fordim (k) {
                if (k > 0) {
                    os << ",";
                }
                os << attr.g[i][k];
            }
            os << " ";
        }
        os.precision(p);
        return os;
    }

    
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

    IntensityForce::IntensityForce(double coeff) : m_Coeff(coeff), m_WorkAtWarpedSpace(false) {

    }

    IntensityForce::~IntensityForce() {

    }

    void IntensityForce::SetImageContext(ImageContext* context) {
        m_ImageContext = context;
    }

    void IntensityForce::SetWorkAtWarpedSpace(bool check) {
        m_WorkAtWarpedSpace = check;
    }

    void IntensityForce::ComputeAttributes(ParticleSystem* system) {
        ParticleSubjectArray& shapes = system->GetSubjects();
        const ParticleSubject& meanSubject = system->GetMeanSubject();

        const int nSubj = shapes.size();
        const int nPoints = shapes[0].GetNumberOfPoints();
        const int nRadius = 1;
        itkcmds::itkImageIO<DoubleImage> io;

        m_attrs.resize(nSubj, nPoints);
        m_attrsMean.set_size(nPoints, Attr::NATTRS);
        warpedImages.resize(nSubj);
        gradImages.resize(nSubj);

        if (m_WorkAtWarpedSpace) {
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
            }
        } else {
            for (int i = 0; i < nSubj; i++) {
                warpedImages[i] = m_ImageContext->GetDoubleImage(i);
                ParticleSubject& subject = shapes[i];
                for (int j = 0; j < nPoints; j++) {
                    forcopy (subject.m_Particles[j].x, subject.m_Particles[j].y);
                }
            }
        }

        for (int i = 0; i < nSubj; i++) {
            ParticleSubject& subject = shapes[i];

            //            char warpedname[128];
            //            sprintf(warpedname, "warped%d_%03d.nrrd", i, system->currentIteration);
            //            io.WriteImageT(warpedname, warpedImages[i]);

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
                Attr& jAttr = m_attrs(i, j);
                for (int k = 0; k < Attr::NATTRS; k++) {
                    IntIndex idx = iiter.GetIndex(k);
                    jAttr.x[k] = warpedImage[idx];
                    fordim (u) {
                        jAttr.g[k][u] = gradImage[idx][u];
                    }
                }
            }
        }
    }

    void IntensityForce::ComputeIntensityForce(ParticleSystem* system) {
        ParticleSubjectArray& shapes = system->GetSubjects();
        const int nSubj = shapes.size();
        const int nPoints = shapes[0].GetNumberOfPoints();

        ComputeAttributes(system);

        // assume all attributes are computed already
        assert(m_attrs.size1() == nSubj && m_attrs.size2() == nPoints);
        
        // column sum
        #pragma omp parallel for
        for (int i = 0; i < nPoints; i++) {
            for (int j = 0; j < Attr::NATTRS; j++) {
                m_attrsMean[i][j] = 0;
                for (int s = 0; s < nSubj; s++) {
                    double attr = m_attrs(s, i).x[j];
                    m_attrsMean[i][j] += attr;
                }
            }
        }
        m_attrsMean /= nSubj;

        // compute mean differences
        #pragma omp parallel for
        for (int i = 0; i < nSubj; i++) {
            for (int j = 0; j < nPoints; j++) {
                for (int k = 0; k < Attr::NATTRS; k++) {
                    m_attrs(i,j).y[k] = m_attrs(i,j).x[k] - m_attrsMean[j][k];
                }
            }
        }
//        __showmatrix(m_attrs, nSubj, nPoints);

        VNLMatrix eye(nSubj, nSubj);
        eye.set_identity();
        for (int i = 0; i < nPoints; i++) {
            // covariance matrix
            VNLMatrix cov(nSubj, nSubj);
            cov.fill(0);
            for (int j = 0; j < nSubj; j++) {
                Attr& jattr = m_attrs(j,i);
                for (int k = j; k < nSubj; k++) {
                    Attr& kattr = m_attrs(k,i);
                    for (int l = 0; l < Attr::NATTRS; l++) {
                        cov[j][k] += jattr.y[l] * kattr.y[l];
                    }
                    cov[k][j] = cov[j][k];
                }
            }
            cov /= (Attr::NATTRS - 1);
            cov = cov + eye;
            VNLMatrix covInverse = vnl_matrix_inverse<double>(cov);

            for (int j = 0; j < nSubj; j++) {
                Attr& jattr = m_attrs(j,i);
                fordim (k) {
                    jattr.f[k] = 0;
                }
                for (int l = 0; l < Attr::NATTRS; l++) {
                    jattr.x[l] = 0;
                    for (int k = 0; k < nSubj; k++) {
                        Attr& kattr = m_attrs(k,i);
                        jattr.x[l] += kattr.y[l] * covInverse[k][j];
                    }
                    fordim (k) {
                        jattr.f[k] += jattr.x[l] * jattr.g[l][k];
                    }
                }
            }
        }
//        __showcmd(attrs, nSubj, nPoints, cout << m_attrs(_,__).f[0] << "," << m_attrs(_,__).f[1]);


        // compute force at subject space
        for (int i = 0; i < nSubj; i++) {
            if (m_WorkAtWarpedSpace) {
                FieldTransformType::Pointer fieldTransform = shapes[i].m_InverseDeformableTransform;
                for (int j = 0; j < nPoints; j++) {
                    double* forcePtr = m_attrs(i,j).f;
                    double* forceOutPtr = m_attrs(i,j).F;
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
            } else {
                for (int j = 0; j < nPoints; j++) {
                    forcopy (m_attrs(i,j).f, m_attrs(i,j).F);
                }
            }
        }
//        __showcmd(attrs, nSubj, nPoints, cout << " " << m_attrs(_,__).F[0] << "," << m_attrs(_,__).F[1]);


        for (int i = 0; i < nSubj; i++) {
            ParticleSubject& subj = shapes[i];
            for (int j = 0; j < nPoints; j++) {
                Particle& par = subj[j];
                VNLVector ff(__Dim);
                fordim (k) {
                    ff[k] = m_attrs(i,j).F[k];
                }
                ff.normalize();
                par.AddForce(ff.data_block(), m_Coeff);
            }
        }
    }
}
