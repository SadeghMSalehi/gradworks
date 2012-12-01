//
//  myParticleDynamics.h
//  ParticlesGUI
//
//  Created by Joohwi Lee on 11/29/12.
//
//

#ifndef __ParticlesGUI__myParticleDynamics__
#define __ParticlesGUI__myParticleDynamics__

#include <iostream>
#include "vnlCommon.h"
#include "itkOptimizerCommon.h"
#include "myImplicitSurfaceConstraint.h"

class ParticleSystem {
    friend class ParticleSystemObserver;

public:
    ParticleSystem(const int nSubj, const int nParticles);
    ~ParticleSystem() {};

    void SetPositions(OptimizerParametersType* params);
    void GetPositions(OptimizerParametersType* params);
    void SetConstraint(myImplicitSurfaceConstraint* constraint);
    void SetHistoryVector(VNLVectorArray* systemHistory);

    void Integrate();

    // functor for integration
    void operator()(const VNLVector &x, VNLVector& dxdt, const double t);

    // functor for observer
    // moved to observer
    // void operator()(const VNLVector &x, const double t);

private:
    ParticleSystem() : m_nDim(0), m_nSubject(0), m_nParticles(0), m_nParams(0) {};

    
    // compute forces between particles
    void UpdateForce(VNLMatrixRef& pos, VNLMatrixRef& vel, VNLMatrix& force);

    // apply constraint on implicit boundaries
    void UpdateConstraint(const VNLVector &x, VNLMatrixRef& dpdt, VNLMatrixRef& dvdt, VNLMatrix& gForce);

    const int m_nDim;
    const int m_nSubject;
    const int m_nParticles;
    const int m_nParams;
    
    double m_Cutoff;
    double m_Sigma2;
    double m_Viscosity;
    double m_Mu;
    double m_COR;

    int m_TimeStep;

    VNLVector m_Status;
    VNLMatrix m_Force;

    VNLVectorArray* m_StatusHistory;
    myImplicitSurfaceConstraint* m_Constraint;
};

class ParticleSystemObserver {
public:
    ParticleSystemObserver(ParticleSystem* system) : m_System(system) {
    }
    void operator()(const VNLVector &x, const double t);
private:
    ParticleSystem* m_System;
};


#endif /* defined(__ParticlesGUI__myParticleDynamics__) */