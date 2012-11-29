//
//  myParticleDynamics.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 11/29/12.
//
//

#include "myParticleDynamics.h"
#include "boost/numeric/odeint.hpp"

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

#define __dist2(x1,y1,x2,y2) ((x1-y1)*(x1-y1)+(x2-y2)*(x2-y2))

static inline double dist2(VNLMatrix& m, int r, int p1, int p2, int d = 2) {
    return __dist2(m[r][d*p1], m[r][d*p1+1], m[r][d*p2], m[r][d*p2+1]);
}

ParticleSystem::ParticleSystem(const int nSubj, const int nParticles): m_nDim(2), m_nSubject(nSubj), m_nParticles(nParticles), m_nParams(nParticles*m_nDim), m_Pos(0,0,NULL), m_Vel(0,0,NULL) {
    m_Cutoff = 15;
    m_Sigma2 = 3;
    m_Mu = 0;
    m_Force.set_size(m_nSubject, m_nParams);
}

void ParticleSystem::SetPositions(OptimizerParametersType* params) {
    m_Status.set_size(m_nSubject*m_nParams*2);
    m_Status.fill(0);
    params->copy_out(m_Status.data_block());
    m_Pos = VNLMatrixRef(m_nSubject, m_nParams, m_Status.data_block());
    m_Vel = VNLMatrixRef(m_nSubject, m_nParams, m_Status.data_block() + m_nSubject * m_nParams);
}

void ParticleSystem::GetPositions(OptimizerParametersType* params) {
    m_Pos.copy_out(params->data_block());
}

void ParticleSystem::UpdateForce() {
    typedef VNLVec2Ref vec_ref;
    typedef VNLVec2 vec;

    const int nDim = vec::SIZE;

    // clear force
    m_Force.fill(0);

    // compute forces between particles
    VNLVector weights(m_nParticles);
    for (int n = 0; n < m_nSubject; n++) {
        for (int i = 0; i < m_nParticles; i++) {
            // reference data
            vec_ref force(&m_Force[n][nDim*i]);
            vec_ref vel(&m_Vel[n][nDim*i]);
            vec_ref pos(&m_Pos[n][nDim*i]);

            for (int j = 0; j < m_nParticles; j++) {
                if (i == j) {
                    // there's no self interaction
                    weights[j] = 0;
                    continue;
                }
                double dij2 = dist2(m_Pos, 0, i, j);
                if (dij2 > m_Cutoff * m_Cutoff) {
                    continue;
                }
                weights[j] = exp(-dij2/(m_Sigma2));
            }
            double sumForce = weights.sum();
            if (sumForce > 0) {
                weights /= sumForce;
            }

            // actual force update
            vec xixj;
            // update force for neighboring particles
            for (int j = 0; j < m_nParticles; j++) {
                if (i == j || weights[j] == 0) {
                    continue;
                }
                vec_ref posj(&m_Pos[n][nDim*i]);
                vec_ref::sub(pos.data_block(), posj.data_block(), xixj.data_block());
                xixj.normalize();
                force += (weights[j] * xixj);
            }

            // dragging force
            force -= (m_Mu * vel);
        }
    }
}

void ParticleSystem::UpdateConstraint() {
    // currently no constraint
}

void ParticleSystem::operator()(const VNLVector &x, VNLVector& dxdt, const double t) {
    VNLMatrixRef pos(m_nSubject, m_nParams, m_Status.data_block());
    VNLMatrixRef vel(m_nSubject, m_nParams, m_Status.data_block()+m_nSubject*m_nParams);

    VNLMatrixRef dpdt(m_nSubject, m_nParams, dxdt.data_block());
    VNLMatrixRef dvdt(m_nSubject, m_nParams, dxdt.data_block()+m_nSubject*m_nParams);

    // dP/dt = V
    dpdt.copy_in(pos.data_block());

    // dV/dt = F/m
    UpdateForce();
    dvdt.copy_in(m_Force.data_block());
}

void ParticleSystem::operator()(const VNLVector &x, const double t) {

}

void ParticleSystem::Integrate() {
    boost::numeric::odeint::integrate((*this), m_Status, 0., 1., 0.1, (*this));
}