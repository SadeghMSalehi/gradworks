//
//  myParticleDynamics.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 11/29/12.
//
//

#include "myParticleDynamics.h"

#define __dist2(x1,y1,x2,y2) ((x1-y1)*(x1-y1)+(x2-y2)*(x2-y2))

static inline double dist2(VNLMatrix& m, int r, int p1, int p2, int d = 2) {
    return __dist2(m[d*p1], m[d*p1+1], m[d*p2], m[d*p2+1]);
}

ParticleSystem::ParticleSystem(const int nSubj, const int nParticles): m_nDim(2), m_nSubject(nSubj), m_nParticles(nParticles), m_nParams(nParticles*m_nDim) {
    m_Cutoff = 15;
    m_Sigma2 = 3;
    m_Pos.set_size(m_nSubject, m_nParams);
    m_Vel.set_size(m_nSubject, m_nParams);
    m_Force.set_size(m_nSubject, m_nParams);
    m_Pos.fill(0);
    m_Vel.fill(0);
    m_Force.fill(0);
}

void ParticleSystem::SetPositions(OptimizerParametersType* params) {
    m_Pos.copy_in(params->data_block());
}

void ParticleSystem::GetPositions(OptimizerParametersType* params) {
    m_Pos.copy_out(params->data_block());
}

void ParticleSystem::UpdateForce() {
    // clear force
    m_Force.fill(0);
    
    // compute forces between particles
    VNLVector forces(m_nParticles);
    for (int n = 0; n < m_nSubject; n++) {
        for (int i = 0; i < m_nParticles; i++) {
            for (int j = 0; j < m_nParticles; j++) {
                if (i == j) {
                    // there's no self interaction
                    forces[j] = 0;
                    continue;
                }
                double dij2 = dist2(m_Pos, 0, i, j);
                if (dij2 > m_Cutoff * m_Cutoff) {
                    continue;
                }
                forces[j] = exp(-dij2/(m_Sigma2));
            }
            double sumForce = forces.sum();
            if (sumForce > 0) {
                forces /= sumForce;
            }
        }
    }
}

void ParticleSystem::UpdateConstraint() {
    // currently no constraint
}

void ParticleSystem::UpdateDerivative() {
    
}

void ParticleSystem::Integrate() {
    
}