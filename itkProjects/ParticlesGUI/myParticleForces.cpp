//
//  myParticleForces.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 1/26/13.
//
//

#include "myParticleCore.h"
#include "myParticleBSpline.h"

namespace pi {
    void InternalForce::ComputeForce(ParticleShapeArray& shapes) {
        const double mu = 1;
        const int nShapes = shapes.size();
        const int nPoints = shapes[0].m_Particles.size();
        for (int k = 0; k < nShapes; k++) {
            ParticleArray& particles = shapes[k].m_Particles;
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

    void EnsembleForce::ComputeMeanShape(ParticleShapeArray& shapes) {
        const int nShapes = shapes.size();
        m_MeanShape.m_SubjId = -1;
        m_MeanShape.m_nPoints = shapes[0].m_nPoints;
        m_MeanShape.Zero();

        for (int i = 0; i < m_MeanShape.m_nPoints; i++) {
            for (int j = 0; j < nShapes; j++) {
                fordim(k) {
                    m_MeanShape[i].x[k] += shapes[j][i].x[k];
                }
            }
            fordim(k) {
                m_MeanShape[i].x[k] /= nShapes;
            }
        }
    }
    
    void EnsembleForce::ComputeForce(ParticleShapeArray& shapes) {
        if (shapes.size() < 2) {
            return;
        }
        const int nPoints = shapes[0].m_nPoints;
        const int nShapes = shapes.size();
        
        ComputeMeanShape(shapes);
        for (int i = 0; i < shapes.size(); i++) {
            ParticleBSpline transform;
            transform.EstimateTransform(m_MeanShape, shapes[i]);
            FieldTransformType::Pointer fieldTransform = transform.GetTransform();
            shapes[i].TransformX2Y(fieldTransform.GetPointer());
            shapes[i].m_Transform = fieldTransform;
        }
        for (int j = 0; j < nPoints; j++) {
            for (int i = 0; i < shapes.size(); i++) {
                fordim(k) {
                    m_MeanShape[j].y[k] += shapes[i][j].y[k];
                }
            }
            fordim(k) {
                m_MeanShape[j].y[k] /= nShapes;
            }
        }

        for (int i = 0; i < shapes.size(); i++) {
            for (int j = 0; j < nPoints; j++) {
                FieldTransformType::InputPointType xPoint;
                FieldTransformType::JacobianType xJac;
                xJac.set_size(2,2);

                double f[4] = { 0, };
                double *x = shapes[i][j].x;
                double *y = shapes[i][j].y;
                double *my = m_MeanShape[j].y;
                fordim(k) {
                    xPoint[k] = x[k];
                }
                shapes[i].m_Transform->ComputeInverseJacobianWithRespectToPosition(xPoint, xJac);
                fordim(k) {
                    f[k] = xJac[0][k]*(y[0]-my[0]) + xJac[1][k]*(y[1]-my[1]) + xJac[2][k]*(y[2]-my[2]);
                }
                shapes[i][j].SubForce(f, 0.1);
            }
        }
    }

}
