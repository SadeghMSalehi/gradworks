//
//  myParticleDynamics.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 11/29/12.
//
//

#include "myParticleDynamics.h"
#include "boost/numeric/odeint.hpp"
#include "iostream"

using namespace std;

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

ParticleSystem::ParticleSystem(const int nSubj, const int nParticles): m_nDim(2), m_nSubject(nSubj), m_nParticles(nParticles), m_nParams(nParticles*m_nDim) {
    m_Cutoff = 15;
    m_Sigma2 = 3*3;
    m_Mu = 1;
    m_COR = 1;
    m_Force.set_size(m_nSubject, m_nParams);
    m_Constraint = NULL;
    m_StatusHistory = NULL;
}

void ParticleSystem::SetPositions(OptimizerParametersType* params) {
    m_Status.set_size(m_nSubject*m_nParams*2);
    m_Status.fill(0);
    params->copy_out(m_Status.data_block());

    VNLMatrixRef gPos(m_nSubject, m_nParams, m_Status.data_block() + m_nSubject * m_nParams);
    gPos.fill(1);
    
}

void ParticleSystem::GetPositions(OptimizerParametersType* params) {
    VNLMatrixRef pos(m_nSubject, m_nParams, m_Status.data_block());
    pos.copy_out(params->data_block());
}

void ParticleSystem::SetHistoryVector(VNLVectorArray* statusHistory) {
    m_StatusHistory = statusHistory;
}

void ParticleSystem::SetConstraint(myImplicitSurfaceConstraint* constraint) {
    m_Constraint = constraint;
}

void ParticleSystem::UpdateForce(VNLMatrixRef& gPos, VNLMatrixRef& gVel, VNLMatrix& gForce) {
    const int nDim = 2;//vec::SIZE;

    // compute forces between particles
    VNLVector weights(m_nParticles);
    for (int n = 0; n < m_nSubject; n++) {
        for (int i = 0; i < m_nParticles; i++) {
            // reference data
            VNLVectorRef force(nDim, &gForce[n][nDim*i]);
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
                weights[j] = exp(-dij*dij/(m_Sigma2));
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
                VNLVectorRef posj(nDim, &gPos[n][nDim*j]);
                VNLCVector::subtract(pos.data_block(), posj.data_block(), xixj.data_block(), nDim);
                xixj.normalize();
                force += (weights[j] * xixj);
            }

            // dragging force
            force -= (m_Mu * vel);

        }
    }
}

void ParticleSystem::UpdateConstraint(const VNLVector &x, VNLMatrixRef& dpdt, VNLMatrixRef& dvdt, VNLMatrix& gForce) {
    
    // currently no constraint
    // boundary constraint
    const int nDim = 2;
    VNLMatrixRef gPos(m_nSubject, m_nParams, (double*) x.data_block());
    VNLMatrixRef gVel(m_nSubject, m_nParams, (double*) x.data_block()+m_nSubject*m_nParams);

    for (int n = 0; n < m_nSubject; n++) {
        for (int i = 0; i < m_nParticles; i++) {
            VNLVectorRef posi(nDim, &gPos[n][nDim*i]);
            VNLVectorRef force(nDim, &gForce[n][nDim*i]);
            VNLVectorRef vel(nDim, &gVel[n][nDim*i]);
            
            if (m_Constraint != NULL) {
                SliceInterpolatorType::ContinuousIndexType idx;
                SliceInterpolatorType::IndexType nidx;
                nidx[0] = idx[0] = posi[0];
                nidx[1] = idx[1] = posi[1];


                double dist = m_Constraint->GetDistance(n, idx);
                GradientType g = m_Constraint->GetGradient(n, nidx);
                VNLVectorRef normal(nDim, g.GetDataPointer());

                //cout << "Dist: " << dist << "; Normal: " << g << "; |Normal|: " << normal.two_norm() << endl;
                if (normal.two_norm() > 0.1) {
                    VNLVectorRef dvdti(nDim, &dvdt[n][nDim*i]);
                    normal.normalize();
                    double velMag = vel.two_norm();
                    VNLVector newVel = velMag * normal - vel;
                    newVel *= 10;
                    //cout << "NewVal: " << newVel << endl;
                    cout << "Velocity Expected to Change : " << vel << " => " << (vel + newVel) << endl;
                    dvdti.set(newVel.data_block());
                }

                /*
                if (dist > -3) {
                    if (normal.two_norm() > 1e-2) {
                        normal.normalize();
                        double forceNormalComponent = dot_product(force, normal);
                        if (forceNormalComponent < 0) {
                            // remove normal component from force
                            force += (forceNormalComponent * normal);
                        }
                        // compute impulsion
                        double velocityNormalComponent = dot_product(vel, normal);
                        if (velocityNormalComponent < 0) {
                            //force += (vel + 2*velocityNormalComponent*normal);
                        }
                    }
                    if (dist >= 0) {
                        normal.normalize();
                        double velocityNormalComponent = dot_product(vel, normal);
                        if (velocityNormalComponent < 0) {
                            //force += (vel + 2*velocityNormalComponent*normal);
                        }
                        for (int k = 0; k < nDim; k++) {
                            force[k] = -vel[k];
                        }
                    }
                }
                 */
            }
        }
    }

}

void ParticleSystem::operator()(const VNLVector &x, VNLVector& dxdt, const double t) {
    VNLMatrixRef pos(m_nSubject, m_nParams, (double*) x.data_block());
    VNLMatrixRef vel(m_nSubject, m_nParams, (double*) x.data_block()+m_nSubject*m_nParams);

    VNLMatrixRef dpdt(m_nSubject, m_nParams, dxdt.data_block());
    VNLMatrixRef dvdt(m_nSubject, m_nParams, dxdt.data_block()+m_nSubject*m_nParams);

    // cout << "Velocity: " << vel << endl;
    // dP/dt = V
    dpdt.copy_in(vel.data_block());

    // update force at time t
    m_Force.fill(0);

    // UpdateForce(pos, vel, m_Force);

    // dV/dt = F/m
    dvdt.copy_in(m_Force.data_block());

    UpdateConstraint(x, dpdt, dvdt, m_Force);
}

void ParticleSystemObserver::operator()(const VNLVector &x, const double t) {
    VNLMatrixRef pos(m_System->m_nSubject, m_System->m_nParams, (double*) x.data_block());
    VNLMatrixRef vel(m_System->m_nSubject, m_System->m_nParams, (double*) x.data_block()+m_System->m_nSubject*m_System->m_nParams);
    if (m_System->m_StatusHistory != NULL) {
        m_System->m_StatusHistory->push_back(x);
    }
    cout << "Time: " << t << endl;
}

void ParticleSystem::Integrate() {
    ParticleSystemObserver observer(this);

    VNLVector status(m_Status);
    // use constant time step
    const int RK4 = 2;
    const int EULER = 1;
    const int RKF45 = 0;

    double t0 = 0;
    double t1 = 100;
    
    int odeMethod = EULER;
    boost::numeric::odeint::euler<VNLVector> eulerStepper;
    boost::numeric::odeint::runge_kutta4<VNLVector> rk4Stepper;
    switch (odeMethod) {
        case RKF45:
            boost::numeric::odeint::integrate((*this), m_Status, t0, t1, .1, observer);
            break;
        case EULER:
            boost::numeric::odeint::integrate_const(eulerStepper, (*this), m_Status, t0, t1, .1, observer);
            break;
        case RK4:
            boost::numeric::odeint::integrate_const(rk4Stepper, (*this), m_Status, t0, t1, .1, observer);
            break;
        default:
            break;
    }

    cout << "History size: " << m_StatusHistory->size() << endl;

}