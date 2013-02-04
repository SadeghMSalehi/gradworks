//
//  myParticleForces.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 1/26/13.
//
//

#include "myParticleCore.h"
#include "myParticleBSpline.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkConstNeighborhoodIterator.h"
#include "armadillo"

namespace pi {
    typedef itk::GradientRecursiveGaussianImageFilter<DoubleImage, VectorImage> GradientFilterType;
    typedef GradientFilterType::OutputPixelType GradientPixelType;
    typedef itk::VectorLinearInterpolateImageFunction<VectorImage> GradientInterpolatorType;
    typedef itk::ConstNeighborhoodIterator<DoubleImage> DoubleImageNeighborhoodIteratorType;
    typedef itk::ConstNeighborhoodIterator<VectorImage> VectorImageNeighborhoodIteratorType;

//    static void ExtractAttributes(DoubleImage::Pointer image, VectorImage::Pointer grad, Particle& par);

    void InternalForce::ComputeForce(ParticleSubjectArray& shapes) {
        const double mu = 1;
        const int nSubjects = shapes.size();
        for (int k = 0; k < nSubjects; k++) {
            const int nPoints = shapes[k].m_Particles.size();
            ParticleArray& particles = shapes[k].m_Particles;
            #pragma omp parallel for
            for (int i = 0; i < nPoints; i++) {
                Particle& pi = particles[i];
                for (int j = i+1; j < nPoints; j++) {
                    Particle& pj = particles[j];
                    ComputeForce(pi, pj);
                }
                double f[__Dim] = { 0, };
                fordim(d) {
                    f[d] = -mu * pi.v[d];
                }
                pi.AddForce(f);
            }
        }
    }
    
    void InternalForce::ComputeForce(Particle &pi, Particle &pj) {
        const double sigma = 15 * 5;
        const double coeff = M_PI_2 / sigma;
        
        double fi[__Dim] = { 0 }, fj[__Dim] = { 0 };
        double dx[__Dim] = { 0 };

        double rij2 = 0;
        pi.Sub(pj, dx);
        fordim(k) {
            rij2 += (dx[k]*dx[k]);
        }
        const double rij = std::sqrt(rij2);

        if (rij <= sigma) {
            fordim(k) {
                dx[k] /= rij;
            }
            const double crij = rij * coeff;
            const double sin1crij = std::sin(crij);
            const double sin2crij = sin1crij * sin1crij;
            fordim(k) {
                fj[k] = fi[k] = dx[k] * (coeff * (1 - (1 / sin2crij)));
            }
            pi.SubForce(fi);
            pj.AddForce(fj);
        }
    }
    

    EnsembleForce::EnsembleForce() {
        
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
    
    void EnsembleForce::ComputeEnsembleForce(ParticleSubjectArray& shapes) {
        if (shapes.size() < 2) {
            return;
        }
        const int nPoints = shapes[0].GetNumberOfPoints();
        const int nSubjects = shapes.size();

        // Compute xMeanShape
        ComputeMeanShape(shapes);
        for (int i = 0; i < shapes.size(); i++) {
            ParticleBSpline transform;
            transform.SetReferenceImage(m_ImageContext->GetLabel(i));
            transform.EstimateTransform(shapes[i], m_MeanShape);
//            transform.EstimateTransform(m_MeanShape, shapes[i]);
            FieldTransformType::Pointer fieldTransform = transform.GetTransform();
            shapes[i].TransformX2Y(fieldTransform.GetPointer());
            shapes[i].m_Transform = fieldTransform;
        }

        // Compute yMeanShape
        fordim(k) {
            for (int j = 0; j < nPoints; j++) {
                for (int i = 0; i < nSubjects; i++) {
                    m_MeanShape[j].y[k] += shapes[i][j].y[k];
                }
                m_MeanShape[j].y[k] /= nSubjects;
            }
        }

        for (int i = 0; i < nSubjects; i++) {
            #pragma omp parallel for
            for (int j = 0; j < nPoints; j++) {
                FieldTransformType::InputPointType xPoint;
                FieldTransformType::JacobianType xJac;
                xJac.set_size(__Dim,__Dim);

                double f[__Dim] = { 0, };
                double *x = shapes[i][j].x;
                double *y = shapes[i][j].y;
                double *my = m_MeanShape[j].y;
                fordim(k) {
                    xPoint[k] = x[k];
                }
                shapes[i].m_Transform->ComputeInverseJacobianWithRespectToPosition(xPoint, xJac);
//                shapes[i].m_Transform->ComputeJacobianWithRespectToPosition(xPoint, xJac);
                fordim(k) {
                    if (__Dim == 3) {
                        f[k] = xJac[0][k]*(y[0]-my[0]) + xJac[1][k]*(y[1]-my[1]) + xJac[2][k]*(y[2]-my[2]);
                    } else if (__Dim == 2) {
                        f[k] = xJac[0][k]*(y[0]-my[0]) + xJac[1][k]*(y[1]-my[1]);
                    }
                }
                shapes[i][j].SubForce(f, 1);
            }
        }
    }

    void EnsembleForce::ComputeImageForce(ParticleSubjectArray &shapes) {
        const int nSubj = shapes.size();
        const int nPoints = shapes[0].GetNumberOfPoints();
        const int nRadius = 2;
        const int nAttrsPerPoint = ::pow(2*nRadius+1, __Dim);
        const int nAttrs = nPoints * nAttrsPerPoint;

        DoubleImageVector warpedImages(nSubj);
        VectorImageVector gradImages(nSubj);

        double* attrs = new double[nSubj * nAttrs];
        double* gradAttrs = new double[nSubj * nAttrs * __Dim];
        double* force = new double[nSubj * nPoints * __Dim];

        // first create a warped image into the mean transform space
        // second create a gradient vector image per subject
        // third extract attributes (features)
        for (int i = 0; i < nSubj; i++) {
            ParticleSubject& subject = shapes[i];
            ParticleBSpline transform;
            transform.SetReferenceImage(m_ImageContext->GetLabel(i));
            transform.EstimateTransform(m_MeanShape, subject);
            FieldTransformType::Pointer fieldTransform = transform.GetTransform();
            subject.m_InverseTransform = fieldTransform;
            warpedImages[i] = transform.WarpImage(m_ImageContext->GetDoubleImage(i));
            GradientFilterType::Pointer grad = GradientFilterType::New();
            grad->SetInput(warpedImages[i]);
            grad->Update();
            gradImages[i] = grad->GetOutput();

            // extract attributes
            DoubleImage::SizeType radius;
            radius.Fill(nRadius);
            DoubleImageNeighborhoodIteratorType iiter(radius, warpedImages[i], warpedImages[i]->GetBufferedRegion());
            VectorImageNeighborhoodIteratorType giter(radius, gradImages[i], gradImages[i]->GetBufferedRegion());

            DoubleImage::IndexType idx;
            #pragma omp parallel for
            for (int j = 0; j < nPoints; j++) {
                Particle& par = subject.m_Particles[j];
                fordim(k) {
                    idx[k] = round(par.y[k]);
                }
                iiter.SetLocation(idx);
                giter.SetLocation(idx);

                double* jAttrs = &attrs[i * nAttrs + j * nAttrsPerPoint];
                double* jAttrsGrad = &gradAttrs[(i * nAttrs + j * nAttrsPerPoint) * __Dim];
                const int numAttrs = iiter.Size();
                assert(numAttrs == nAttrsPerPoint);
                for (int k = 0; k < nAttrsPerPoint; k++) {
                    double pixel = iiter.GetPixel(k);
                    VectorType grad = giter.GetPixel(k);
                    jAttrs[k] = pixel;
                    for (int m = 0; m < __Dim; m++) {
                        jAttrsGrad[m] = grad[m];
                    }
                }
            }
        }


        // column sum
        double* sumAttrs = new double[nPoints * nAttrsPerPoint];
        #pragma omp parallel for
        for (int j = 0; j < nAttrs; j++) {
            sumAttrs[j] = 0;
            for (int i = 0; i < nSubj; i++) {
                sumAttrs[j] += attrs[i * nAttrs + j];
            }
        }

        // compute mean differences
        for (int i = 0; i < nSubj; i++) {
            #pragma omp parallel for
            for (int j = 0; j < nAttrs; j++) {
                attrs[i * nAttrs + j] -= (sumAttrs[j] / nSubj);
            }
        }

        // compute force direction
        memset(force, 0, nSubj * nPoints * __Dim * sizeof(double));
        for (int i = 0; i < nSubj; i++) {
            #pragma omp parallel for
            for (int j = 0; j < nPoints; j++) {
                double* forcePtr = &force[(i * nPoints + j) * __Dim];
                for (int k = 0; k < nAttrsPerPoint; k++) {
                    double attr = attrs[i * nAttrs + j * nAttrsPerPoint];
                    double* gradAttrPtr = &gradAttrs[(i * nAttrs + j * nAttrsPerPoint + k) * __Dim];
                    for (int m = 0; m < __Dim; m++) {
                        forcePtr[m] += attr * gradAttrPtr[m];
                    }
                }
            }
        }


        // compute force at subject space
        for (int i = 0; i < nSubj; i++) {
            FieldTransformType::Pointer transform = shapes[i].m_InverseTransform;
            for (int j = 0; j < nPoints; j++) {
                double* forcePtr = &force[(i * nPoints + j) * __Dim];
                FieldTransformType::InputPointType x;
                fordim(k) {
                    x[k] = shapes[i][j].x[k];
                }
                FieldTransformType::JacobianType jac;
                jac.set_size(__Dim, __Dim);
                transform->ComputeInverseJacobianWithRespectToPosition(x, jac);
                fordim(k) {
                    double ff = 0;
                    if (__Dim == 3) {
                        ff = jac[0][k]*forcePtr[k] + jac[1][k]*forcePtr[k] + jac[2][k]*forcePtr[k];
                    } else if (__Dim == 2) {
                        ff = jac[0][k]*forcePtr[k] + jac[1][k]*forcePtr[k];
                    }
                    forcePtr[k] = ff;
                }
            }
        }

        for (int i = 0; i < nSubj; i++) {
            ParticleSubject& subj = shapes[i];
            for (int j = 0; j < nPoints; j++) {
                Particle& par = subj[j];
                VNLVector ff(__Dim);
                fordim (k) {
                    const int forceIdx = (i * nPoints + j) * __Dim;
                    ff[k] = force[forceIdx + k];
                }
                ff.normalize();
                par.SubForce(ff.data_block(), 1);
            }
        }


        delete[] force;
        delete[] attrs;
        delete[] gradAttrs;
    }
}
