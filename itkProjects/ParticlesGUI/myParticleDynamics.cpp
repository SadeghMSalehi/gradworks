//
//  myParticleDynamics.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 11/29/12.
//
//

#include <cmath>
#include "myParticleDynamics.h"
#include "boost/numeric/odeint.hpp"
#include "vnl/algo/vnl_symmetric_eigensystem.h"
#include "iostream"
#include "myImageParticlesAlgorithm.h"
#include "itkThinPlateSplineKernelTransform.h"
#include "myImageTransform.h"
#include <QElapsedTimer>

#define __dist2(x1,y1,x2,y2) ((x1-y1)*(x1-y1)+(x2-y2)*(x2-y2))

// VNL-Boost compatibility functions
namespace boost {
    namespace numeric {
        namespace odeint {
            typedef VNLVector state_type;
            
            template <>
            struct is_resizeable<state_type> {
                typedef boost::true_type type;
                static const bool value = type::value;
            };
            
            template<>
            struct same_size_impl<state_type, state_type> { // define how to check size
                static bool same_size(const state_type &v1, const state_type &v2) {
                    return v1.size() == v2.size();
                }
            };
            
            template<>
            struct resize_impl<state_type , state_type>
            { // define how to resize
                static void resize( state_type &v1 ,
                                   const state_type &v2) {
                    v1.set_size(v2.size());
                }
            };
        }
    }
}

namespace my {
    
    using namespace std;
    const static int nDim = 2;
    
    typedef itk::ThinPlateSplineKernelTransform<double,nDim> TPSTransformType;
    
    static inline double dist2(VNLMatrix& m, int r, int p1, int p2, int d = 2) {
        return __dist2(m[r][d*p1], m[r][d*p1+1], m[r][d*p2], m[r][d*p2+1]);
    }
    
    ParticleSystem::ParticleSystem(const int nSubj, const int nParticles):
    m_nDim(2),
    m_nSubjects(nSubj),
    m_nParticles(nParticles),
    m_nParams(nParticles*m_nDim),
    m_Pos(0,0,NULL),
    m_Vel(0,0,NULL),
    m_dpdt(0,0,NULL),
    m_dvdt(0,0,NULL) {
        m_Cutoff = 15;
        m_Sigma = 15;
        m_Mu = 1;
        m_COR = 1 ;
        m_Force.set_size(m_nSubjects, m_nParams);
        m_Constraint = NULL;
        m_StatusHistory = NULL;
        m_Callback = NULL;
        m_Context = NULL;
        m_ForceType = 0;
    }
    
    void ParticleSystem::SetPositions(OptimizerParametersType* params) {
        m_Status.set_size(m_nSubjects*m_nParams*2);
        m_Status.fill(0);
        params->copy_out(m_Status.data_block());
        
        //    VNLMatrixRef gPos(m_nSubject, m_nParams, m_Status.data_block() + m_nSubject * m_nParams);
        //    gPos.fill(1);
        
    }
    
    void ParticleSystem::GetPositions(OptimizerParametersType* params) {
        VNLMatrixRef pos(m_nSubjects, m_nParams, m_Status.data_block());
        pos.copy_out(params->data_block());
    }
    
    void ParticleSystem::SetHistoryVector(VNLVectorArray* statusHistory) {
        m_StatusHistory = statusHistory;
    }
    
    void ParticleSystem::SetCostHistoryVector(STDDoubleArray* costHistory) {
        m_CostHistory = costHistory;
    }
    
    void ParticleSystem::SetEventCallback(EventCallback* callback) {
        m_Callback = callback;
    }
    
    void ParticleSystem::SetContext(ImageParticlesAlgorithm* context) {
        m_Context = context;
        
        m_Cutoff = context->GetProperty().GetDouble("cutoffDistance", 15.0);
        m_Sigma  = context->GetProperty().GetDouble("sigma", 3.0);
        m_GradientScale = context->GetProperty().GetDouble("gradientScale", 10.0);

        m_Mu = context->GetProperty().GetDouble("Mu", 1);
        m_COR = context->GetProperty().GetDouble("COR", 0.5);
        m_Options.useSurfaceForce = context->GetProperty().GetBool("actionUseSurfaceForce", true);
        m_Options.useImageForce = context->GetProperty().GetBool("actionUseImageForce", true);
        m_Options.useBoundaryCondition = context->GetProperty().GetBool("actionUseBoundaryConditions", true);
        m_Options.useEnsembleForce = context->GetProperty().GetBool("actionUseEnsembleForce", true);
        m_Options.useAdaptiveSampling = context->GetProperty().GetBool("actionUseAdaptiveControl", false);

        m_Constraint = context->GetConstraint();
    }
    
    void ParticleSystem::UpdateSurfaceForce() {
        VNLMatrixRef& gPos = m_Pos;
        VNLMatrixRef& gVel = m_Vel;
        VNLMatrix& gForce = m_Force;
        
        const int nDim = 2;//vec::SIZE;
        int nSubj = (m_Options.applySurfaceEntropyToFirstOnly ? 1 : m_nSubjects);

        const bool useAdaptiveSampling = m_Options.useAdaptiveSampling;

        // use cotangent-based force
        m_ForceType = 1;

        // compute forces between particles
        VNLVector weights(m_nParticles);
        for (int n = 0; n < nSubj; n++) {
            SliceInterpolatorType::Pointer kappaIntp;
            if (useAdaptiveSampling) {
                kappaIntp = m_Context->GetAttributeInterpolators()->at(n);
            }
            for (int i = 0; i < m_nParticles; i++) {
                // reference data
                VNLVectorRef iForce(nDim, &gForce[n][nDim*i]);
                VNLVectorRef iVel(nDim, &gVel[n][nDim*i]);
                VNLVectorRef iPos(nDim, &gPos[n][nDim*i]);
                
                ContinuousIndexType iIdx;
                iIdx[0] = iPos[0];
                iIdx[1] = iPos[1];

                // compute energy and force
                if (m_ForceType == 0) {
                    const double sigma2 = m_Sigma * m_Sigma;

                    // iteration over particles
                    // may reduce use symmetric properties
                    for (int j = 0; j < m_nParticles; j++) {
                        if (i == j) {
                            // there's no self interaction
                            weights[j] = 0;
                        } else {
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
                            double dij = (iPos-jPos).two_norm();
                            if (dij > m_Cutoff) {
                                weights[j] = 0;
                            } else {
                                weights[j] = exp(-dij*dij*kappa/(sigma2));
                                // debug: check kappa is different between neighbors
                                //                        if (i == 10) {
                                //                            cout << "Distance: " << dij << "; Kappa: " << kappa << endl;
                                //                        }
                            }
                        }
                    }
                    double sumForce = weights.sum();
                    if (sumForce > 0) {
                        weights /= sumForce;
                    }

                    // actual force update
                    VNLVec2 xixj;
                    // update force for neighboring particles
                    for (int j = 0; j < m_nParticles; j++) {
                        if (i == j || weights[j] == 0) {
                            continue;
                        }
                        VNLVectorRef jPos(nDim, &gPos[n][nDim*j]);
                        xixj.normalize();
                        iForce += (weights[j] * xixj);
                    }
                } else if (m_ForceType == 1) {
                    const double sigma = m_Sigma * 5;
                    const double coeff = M_PI_2 / sigma;
                    double energy = 0;
                    VNLVec2 dEdr;
                    for (int j = 0; j < m_nParticles; j++) {
                        VNLVectorRef jPos(nDim, &gPos[n][nDim*j]);
                        // energy derivative
                        VNLCVector::subtract(iPos.data_block(), jPos.data_block(), dEdr.data_block(), nDim);
                        double rij = dEdr.two_norm();
                        dEdr.normalize();
                        if (rij > sigma) {
                            // no contribution
                            continue;
                        } else {
                            rij *= coeff;
                            if (rij > 0) {
                                energy = 1 / std::tan(rij) + rij - M_PI_2;
                                double sin2rij = std::sin(rij);
                                sin2rij *= sin2rij;
                                dEdr *= (coeff * (1 - 1 / sin2rij));
                                iForce -= dEdr;
                            }
                        }
                    }
                }
            }
        }
    }
    
    
    
    void ParticleSystem::UpdateGravityForce() {
        VNLMatrix& gForce = m_Force;
        const int nDim = 2;
        for (int n = 0; n < m_nSubjects; n++) {
            for (int i = 0; i < m_nParticles; i++) {
                gForce[n][nDim*i+1] = 0.98;
            }
        }
    }
    
    void ParticleSystem::EstimateRigidTransform(VNLMatrixRef& gPos, VNLMatrixArray& transforms, VNLMatrixArray& jacobians) {
        if (m_nParticles < 5) {
            cout << "Not enough particles for procrustes estimation" << endl;
            return;
        }
        
        VNLVector gPosMean(m_nParams);
        vnl_row_mean(gPos, gPosMean);
        
        // move target points to origin center
        VNLMatrix targetPoints(gPosMean.data_block(), m_nParticles, nDim);
        VNLVec2 targetCentroid;
        vnl_row_mean(targetPoints, targetCentroid);
        vnl_row_subtract(targetPoints, targetCentroid, targetPoints);
        
        for (int n = 0; n < gPos.rows(); n++) {
            // move source points to origin center
            VNLMatrix sourcePoints(gPos[n], m_nParticles, nDim);
            VNLVec2 sourceCentroid;
            vnl_row_mean(sourcePoints, sourceCentroid);
            vnl_row_subtract(sourcePoints, sourceCentroid, sourcePoints);
            
            
            double rotationParam;
            VNLVec2 translationParam = targetCentroid - sourceCentroid;
            VNLMatrix cov = sourcePoints.transpose() * targetPoints;
            
            // debug: if the covariance is computed correctly
            //    std::cout << "COV: " << cov << std::endl;
            vnl_svd<double> svd(cov);
            VNLMatrix rotationMatrix = svd.V() * svd.U().transpose();
            
            //    std::cout << rotationMatrix << std::endl;
            
            // debug: atan2 is more stable than acos
            rotationParam = atan2(rotationMatrix[1][0], rotationMatrix[0][0]);
            //    std::cout << rotationParam << std::endl;
            //
            //    std::cout << "acos(1) = " << acos(1) << std::endl;
            
            VNLVector offset = rotationMatrix * translationParam;
            
            VNLMatrix transformMatrix(2,3);
            transformMatrix[0][2] = offset[0];
            transformMatrix[1][2] = offset[1];
            transformMatrix.update(rotationMatrix);
            
            VNLMatrix jacobian(2,2);
            jacobian.update(rotationMatrix);
            
            //        vnl_identity(transformMatrix);
            //        vnl_identity(jacobian);
            
            transforms.push_back(transformMatrix);
            jacobians.push_back(jacobian);
            
            cout << "Estimated Rigid Transform: " << transformMatrix << endl;
        }
    }
    
    void ParticleSystem::ApplyMatrixOperation(const double* posIn, const VNLMatrix& matrix, double* posOut) {
        double *pIn = (double*) posIn;
        double *pOut = posOut;
        for (int i = 0;i < m_nParticles; i++) {
            VNLVec3 posi;
            posi.copy_in(pIn);
            posi[2] = 1;
            VNLVectorRef tposi(nDim, pOut);
            vnl_matrix_x_vector(matrix, posi, tposi);
            pIn += nDim;
            pOut += nDim;
        }
    }
    
    void ParticleSystem::UpdateEnsembleForce()
    {
        VNLMatrixRef& gPos = m_Pos;
        
        // estimate transform from subjN to subj1
        VNLMatrixArray transforms;
        VNLMatrixArray jacobians;
        //    cout << "Position before estimation: " << gPos << endl;
        EstimateRigidTransform(gPos, transforms, jacobians);
        //    cout << "Position after estimation: " << gPos << endl;
        // transform particles onto subj1 space
        VNLMatrix tPos(gPos);
        for (int n = 0; n < m_nSubjects; n++) {
            ApplyMatrixOperation(gPos[n], transforms[n], tPos[n]);
        }
        
        // debug: check if transform is correct
        //    cout << "Before: " << gPos << endl;
        //    cout << "After: " << tPos << endl;
        
        // aggregate attribute data and compute mean
        VNLMatrix data(tPos.rows(), m_nParams);
        data.update(tPos);
        
        VNLVector dataMean(data.cols());
        vnl_row_mean(data, dataMean);
        
        //    cout << "Data Mean: " << dataMean << endl;
        
        // move data to center
        for (int i = 0; i < m_nSubjects; i++) {
            VNLCVector::subtract(data[i], dataMean.begin(), data[i], data.cols());
        }
        
        // debug: positional entropy should work without image params
        //    cout << "nPaarms: " << nParams << endl;
        //    cout << "Data: " << data << endl;
        VNLMatrix cov = data * data.transpose();
        
        // relaxation parameter for singular matrix
        // this produce pseudo-inverse matrix,
        // if alpha is zero, the inverse matrix will have really high numbers
        double alpha = 1;
        for (int i = 0; i < cov.rows(); i++) {
            cov[i][i] += alpha;
        }
        
        // debug: check covariance matrix
        cout << "Data COV: " << cov << endl;
        
        // cost function is the sum of log of eigenvalues
        vnl_symmetric_eigensystem<double> eigen(cov);
        
        // debug: eigenvalues for singular matrix; still producing eigenvalue with zero
        //    cout << "Eigenvalues: " << eigen.D << endl;
        //    cout << "Eigenvectors: " << eigen.V << endl;
        double cost = 0;
        for (int i = 0; i < eigen.D.size(); i++) {
            if (eigen.D[i] > 1) {
                cost += log(eigen.D[i]);
            }
        }
        
        // debug: let's see difference between image gradient and position gradient
        //    cout << "Covariance:" << cov << endl;
        VNLMatrix covInverse = vnl_matrix_inverse<double>(cov);
        VNLMatrix grad = covInverse * data;
        // debug: gradient is the direction minimizing covariance matrix
        //    cout << "Ensemble COV:" << cov << endl;
        //    cout << "Ensemble InverseCOV: " << covInverse << endl;
        //    cout << "gradient: " << grad << endl;
        if (grad.has_nans()) {
            cout << "Gradient has Nans: " << grad << endl;
        }
        
        // multiply jacobian of the function to gradient
        for (int n = 0; n < m_nSubjects; n++) {
            ApplyMatrixOperation(grad[n], jacobians[n], grad[n]);
            grad *= -10;
            //        cout << "Applying gradient: " << grad << endl;
            VNLCVector::add(gPos[n], grad[n], gPos[n], m_nParams);
        }
    }
    
    void ParticleSystem::UpdateImageForce() {
        int nTmpSubjects = 2;
        ImageParticlesAlgorithm::InterpolatorList* interpolators = m_Context->GetImageInterpolators();
        ImageParticlesAlgorithm::GradientInterpolatorList* gradientInterpolators = m_Context->GetGradientInterpolators();
        
        VNLMatrixRef& gPos = m_Pos;
        VNLMatrix wPos(nTmpSubjects, m_nParams);
        for (int n = 0; n < nTmpSubjects; n++) {
            ImageContainer::Pointer image = m_Context->GetImage(n);
            image->TransformToPhysicalPoints(m_nParticles, gPos[n], wPos[n]);
        }
        
        VNLMatrix imageAttributes(nTmpSubjects, m_nParticles);
        for (int n = 0; n < nTmpSubjects; n++) {
            // iterate over particles and sample image intensity
            SliceInterpolatorType::Pointer interp = interpolators->at(n);
            for (int i = 0; i < m_nParticles; i++) {
                VNLVectorRef iPos(nDim, &m_Pos[n][nDim*i]);
                VNLVectorRef iVel(nDim, &m_Vel[n][nDim*i]);
                VNLVectorRef iForce(nDim, &m_Force[n][nDim*i]);
                
                ContinuousIndexType iIdx;
                copyArray(iIdx, iPos);
                
                imageAttributes[n][i] = interpolators->at(n)->EvaluateAtContinuousIndex(iIdx);
            }
        }
        
        // debug: check if intensity sampled correctly
        //    cout << "Image Data: " << endl;
        //    cout << imageAttributes << endl;
        
        // compute covariance matrix and gradient for minimization
        // jacobian requires image intensity gradient with respect to xy coordinate
        VNLVector meanAttributes(m_nParticles);
        VNLMatrix normalizedAttributes(nTmpSubjects, m_nParticles);
        
        vnl_row_mean(imageAttributes, meanAttributes);
        vnl_row_subtract(imageAttributes, meanAttributes, normalizedAttributes);
        
        // debug: check if covariance is correct
        //    cout << "Normalized Attrs: " << normalizedAttributes << endl;
        
        VNLMatrix cov = normalizedAttributes * normalizedAttributes.transpose();
        VNLMatrix covIdentity(cov);
        covIdentity.set_identity();
        cov = cov + covIdentity;
        
        VNLMatrix covInv = vnl_matrix_inverse<double>(cov);
        VNLMatrix grad = covInv * normalizedAttributes;
        
        // debug: check if covariance is correct
        //    cout << "CoV: " << cov << endl;
        //    cout << "Inverse of CoV: " << covInv << endl;
        //    cout << "Gradient: " << grad << endl;
        
        // jacobian * grad
        VNLMatrixRef gForce(nTmpSubjects, m_nParams, m_Force[0]);
        for (int n = 0; n < nTmpSubjects; n++) {
            int kParticle = 0;
            for (int i = 0; i < m_nParams; i+=2) {
                ContinuousIndexType iIdx;
                iIdx[0] = gPos[n][i];
                iIdx[1] = gPos[n][i+1];
                GradientType g = gradientInterpolators->at(n)->EvaluateAtContinuousIndex(iIdx);
                gForce[n][i] -= m_GradientScale * g[0] * grad[n][kParticle];
                gForce[n][i+1] -= m_GradientScale * g[1] * grad[n][kParticle];
                kParticle ++;
            }
        }
    }
    
    void ParticleSystem::UpdateBSplineEnsemble() {
        VNLMatrixRef& xPos = m_Pos;

        if (xPos.has_nans()) {
            cout << "NaNs: " << xPos << endl;
        }
        // ignore conversion between index and physical space
        VNLVector xMeanPos(m_nParams);
//        vnl_row_mean(xPos, xMeanPos);

        // debug: what if we use the first subject as a common space?
        xMeanPos.copy_in(xPos[0]);

        // transform from subject's space to mean space
        std::vector<FieldTransformType::Pointer> bsplineTransformArray;
        for (int n = 0; n < m_nSubjects; n++) {
            SliceType::Pointer refImage = m_Context->GetImage(n)->GetSlice();
            my::BSplineRegistration bReg;
            bReg.SetLandmarks(m_nParticles, xPos[n], xMeanPos.data_block());
            bReg.SetReferenceImage(refImage);
            bReg.SetNumberOfControlPoints(16);
            bReg.Update();
            FieldTransformType::Pointer transformToMean = bReg.GetTransform();
            bsplineTransformArray.push_back(transformToMean);
        }
        
        // compute particles' position at the common space
        //
        VNLMatrix tPos(m_nSubjects, m_nParams);
        tPos.fill(0);

        VNLMatrix jacobPos(m_nSubjects, m_nParams);
        FieldTransformType::InputPointType inPoint;
        FieldTransformType::OutputPointType outPoint;
        for (int n = 0; n < m_nSubjects; n++) {
            if (bsplineTransformArray[n].IsNull()) {
                // if bspline is not computable,
                // use identity transform
                for (int i = 0; i < m_nParticles; i++) {
                    tPos[n][2*i] = xPos[n][2*i];
                    tPos[n][2*i+1] = xPos[n][2*i+1];
                }
            } else {
                for (int i = 0; i < m_nParticles; i++) {
                    inPoint[0] = xPos[n][2*i];
                    inPoint[1] = xPos[n][2*i+1];
                    outPoint = bsplineTransformArray[n]->TransformPoint(inPoint);
                    tPos[n][2*i] = outPoint[0];
                    tPos[n][2*i+1] = outPoint[1];
                }
            }
        }
        
        // gPos: Xi, xPos: T(Xi)
        // now tPos must move towards mean point
        VNLMatrix yPos(tPos);
        vnl_row_center(yPos);

        FieldTransformType::InputPointType xPoint;
        FieldTransformType::JacobianType xJac;
        
        VNLMatrixRef gForce(m_nSubjects, m_nParams, m_Force[0]);
        for (int n = 0; n < m_nSubjects; n++) {
            if (bsplineTransformArray[n].IsNull()) {
                continue;
            }

            // ignore covariance factor
            // because the minimum is same,
            // but it is better to use distance function given in
            // Miriah's paper
            for (int i = 0; i < m_nParticles; i++) {
                xPoint[0] = xPos[n][2*i];
                xPoint[1] = xPos[n][2*i+1];
                bsplineTransformArray[n]->ComputeJacobianWithRespectToPosition(xPoint, xJac);
                double dGdX = xJac[0][0] * yPos[n][2*i] + xJac[1][0] * yPos[n][2*i+1];
                double dGdY = xJac[0][1] * yPos[n][2*i] + xJac[1][1] * yPos[n][2*i+1];
                gForce[n][2*i] -= dGdX;
                gForce[n][2*i+1] -= dGdY;
            }
        }
        
        /*
         // compute covariance and its inverse
         VNLMatrix xCov;
         vnl_covariance(xPos, xCov);
         double alpha = 1;
         vnl_add_diagonal(xCov, alpha);
         
         // inverse of covariance
         VNLMatrix xCovInv = vnl_matrix_inverse<double>(xCov);
         
         // compute the gradient of positional entropy with respect to transformed position
         VNLMatrix pGpP = xCovInv * xPos;
         VNLMatrix pGpX(m_nSubjects, m_nParticles);
         
         // compute the gradient of positional entropy with respect to transformed position
         FieldTransformType::InputPointType xPoint;
         FieldTransformType::JacobianType xJac;
         for (int n = 0; n < m_nSubjects; n++) {
         for (int i = 0; i < m_nParticles; i++) {
         xPoint[0] = gPos[n][2*i];
         xPoint[1] = gPos[n][2*i+1];
         bsplineTransformArray[n]->ComputeJacobianWithRespectToPosition(xPoint, xJac);
         }
         }
         */
        
    }
    
    void ParticleSystem::UpdateKernelTransform() {
        VNLMatrixRef& gPos = m_Pos;
        
        VNLMatrix wPos(m_nSubjects, m_nParams);
        for (int n = 0; n < m_nSubjects; n++) {
            m_Context->GetImage(n)->TransformToPhysicalPoints(m_nParticles, gPos[0], wPos[0]);
        }
        
        
        QElapsedTimer timer;
        timer.start();
        
        my::ImageTransform transformer;
        my::KernelTransformPointer transform = transformer.CreateKernelTransform(0, m_nParticles, gPos[0], wPos[0]);
        
        cout << "Transform Time: " << timer.elapsed() << " ms" << endl;
    }
    
    // Apply boundary constraint
    void ParticleSystem::ApplyBoundaryConditions() {
        
        // boundary constraint
        const int nDim = 2;
        
        int nSubj = m_Options.applyBoundaryConditionToFirstOnly ? 1 : m_nSubjects;
        nSubj = m_nSubjects;
        
        // iterate over all subjects
        for (int n = 0; n < nSubj; n++) {
            // iterate over all particles
            for (int i = 0; i < m_nParticles; i++) {
                // input data
                VNLVectorRef posi(nDim, &m_Pos[n][nDim*i]);
                VNLVectorRef forcei(nDim, (double*) &m_Force[n][nDim*i]);
                VNLVectorRef iVel(nDim, &m_Vel[n][nDim*i]);
                
                // output data
                VNLVectorRef dpdti(nDim, &m_dpdt[n][nDim*i]);
                VNLVectorRef dvdti(nDim, &m_dvdt[n][nDim*i]);
                
                // constraint should be available
                if (m_Constraint != NULL) {
                    SliceInterpolatorType::ContinuousIndexType idx;
                    SliceInterpolatorType::IndexType nidx;
                    nidx[0] = idx[0] = posi[0];
                    nidx[1] = idx[1] = posi[1];
                    
                    if (!m_Constraint->IsInsideRegion(n, nidx)) {
                        iVel.fill(0);
                        continue;
                    }
                    
                    // this is to move outside particles to the nearest boundary
                    double dist = m_Constraint->GetDistance(n, idx);
                    if (dist >= 0 && true) {
                        ImplicitSurfaceConstraint::OffsetType offset = m_Constraint->GetOutsideOffset(n, nidx);
                        for (int k = 0; k < nDim; k++) {
                            posi[k] += offset[k];
                        }
                    }
                    
                    // compute on new position
                    nidx[0] = idx[0] = posi[0];
                    nidx[1] = idx[1] = posi[1];
                    
                    //                if (m_Constraint->GetOutsideOffset(n, nidx))
                    
                    // check boundaries
                    GradientType g = m_Constraint->GetGradient(n, idx);
                    VNLVectorRef normal(nDim, g.GetDataPointer());
                    double normalMagnitude = normal.two_norm();
                    
                    // normalize to compute direction
                    //                normal.normalize();
                    
                    // boundary conditions
                    // velocity set to zero
                    // the normal component of the force set to zero
                    //                double velMag = veli.two_norm();
                    
                    // distance and gradient sometimes doesn't coincide.
                    if (normalMagnitude > 0.4) {
                        normal.normalize();
                        
                        // velocity should be zero toward normal direction
                        double normalSpeed = dot_product(iVel, normal);
                        if (normalSpeed < 0) {
                            VNLVector newVelocity = iVel - 2 * normalSpeed * normal;
                            newVelocity *= m_COR;
                            // how to know current timestep?
                            //                        newVelocity /= 0.1;
                            dvdti.fill(0);
                            iVel.copy_in(newVelocity.data_block());
                            //                        cout << "Velocity: " << veli << " => " << newVelocity << "; Speed: " << normalSpeed << "; Normal: " << normal << endl;
                        }
                        
                        // remove normal term of the current force
                        double normalForce = dot_product(normal, forcei);
                        if (normalForce < 0) {
                            VNLVector newForce = forcei + normalForce * normal;
                            dvdti.copy_in(newForce.data_block());
                            //                        cout << "Force: " << forcei << ", " << newForce << "; grad: " << normal << endl;
                        }
                    }
                }
            }
        }
        
    }


    void ParticleSystem::UpdateDraggingForce() {
        VNLMatrixRef iForce(m_nSubjects, m_nParams, m_Force[0]);
        for (int n = 0; n < m_nSubjects; n++) {
            for (int i = 0; i < m_nParams; i += SDim) {
                VNLVectorRef iVel(SDim, &m_Vel[n][i]);
                VNLVectorRef iForce(SDim, &m_Force[n][i]);
                iForce[0] -= m_Mu * iVel[0];
                iForce[1] -= m_Mu * iVel[1];
            }
        }
    }
    
    
    // integration function
    void ParticleSystem::operator()(const VNLVector &x, VNLVector& dxdt, const double t) {
        const int nVelocityOffset = m_nSubjects * m_nParams;
        
        m_Pos.reinit(m_nSubjects, m_nParams, (double*) &x[0]);
        m_Vel.reinit(m_nSubjects, m_nParams, (double*) &x[nVelocityOffset]);
        m_dpdt.reinit(m_nSubjects, m_nParams, (double*) &dxdt[0]);
        m_dvdt.reinit(m_nSubjects, m_nParams, &dxdt[nVelocityOffset]);
        
        m_Force.fill(0);
        
        // update forces at time t
        if (m_Options.useEnsembleForce) {
            UpdateBSplineEnsemble();
        }
        
        if (m_Options.useImageForce) {
            UpdateImageForce();
        }
        
        if (m_Options.useSurfaceForce) {
            UpdateSurfaceForce();
        }

        UpdateDraggingForce();

        if (m_Options.useParticlePhysics) {
            // dP/dt = V
            m_dpdt.copy_in(m_Vel[0]);
            // dV/dt = F/m
            m_dvdt.copy_in(m_Force[0]);
        } else {
            // gradient descent
            m_dpdt.copy_in(m_Force.data_block());
            m_dvdt.fill(0);
        }
        
        if (m_Options.useBoundaryCondition) {
            ApplyBoundaryConditions();
        }
    }
    
    
    // observer function
    void ParticleSystem::operator()(const VNLVector &x, const double t) {
        const int nDim = 2;
        VNLMatrixRef gPos(m_nSubjects, m_nParams, (double*) x.data_block());
        VNLMatrixRef gVel(m_nSubjects, m_nParams, (double*) x.data_block()+m_nSubjects*m_nParams);
        
        double cost = 0;
        // compute forces between particles
        VNLVector weights(m_nParticles);
        const double sigma2 = m_Sigma*m_Sigma;
        for (int n = 0; n < m_nSubjects; n++) {
            for (int i = 0; i < m_nParticles; i++) {
                // reference data
                VNLVectorRef vel(nDim, &gVel[n][nDim*i]);
                VNLVectorRef pos(nDim, &gPos[n][nDim*i]);
                
                for (int j = 0; j < m_nParticles; j++) {
                    if (i == j) {
                        // there's no self interaction
                        weights[j] = 0;
                        continue;
                    }
                    VNLVectorRef posj(2, &gPos[n][nDim*j]);
                    double dij = (pos-posj).two_norm();
                    if (dij > m_Cutoff) {
                        weights[j] = 0;
                        continue;
                    }
                    weights[j] = exp(-dij*dij/(sigma2));
                }
                cost += weights.sum();
            }
        }
        
        double xy[2];
        xy[0] = t;
        xy[1] = cost;
        if (m_Callback != NULL) {
            m_Callback->EventRaised(0xADDCEC, 0, NULL, xy);
        }
        
        
        if (m_StatusHistory != NULL) {
            m_StatusHistory->push_back(x);
            m_CostHistory->push_back(cost);
        }
        cout << "Time: " << t << endl;
    }
    
    void ParticleSystem::Integrate(VNLVector& status, double dt, double t0, double t1, int odeMethod) {
        //    ParticleSystemObserver observer(this);
        // use constant time step
        const int RK4 = 2;
        const int EULER = 1;
        const int RKF45 = 0;

        m_Timer.start();
        
        boost::numeric::odeint::euler<VNLVector> eulerStepper;
        boost::numeric::odeint::runge_kutta4<VNLVector> rk4Stepper;
        switch (odeMethod) {
            case RKF45:
                boost::numeric::odeint::integrate((*this), status, t0, t1, dt, (*this));
                break;
            case EULER:
                boost::numeric::odeint::integrate_const(eulerStepper, (*this), status, t0, t1, dt, (*this));
                break;
            case RK4:
                boost::numeric::odeint::integrate_const(rk4Stepper, (*this), status, t0, t1, dt, (*this));
                break;
            default:
                break;
        }

        qint64 elapsedTime = m_Timer.elapsed();
        cout << "Elapsed Time: " << elapsedTime << " ms (" << (elapsedTime / (t1- t0) * dt) << ") ms" << endl;
    }
    
    void ParticleSystem::Integrate() {
        double dt = m_Context->GetProperty().GetDouble("timeStep", 0.1);
        Integrate(m_Status, dt, 0, m_Context->GetProperty().GetInt("numberOfIterations", 100) * dt);
    }
    
}