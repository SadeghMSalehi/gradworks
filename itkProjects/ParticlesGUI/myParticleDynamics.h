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

class ParticleSystem {
public:
    ParticleSystem(const int nSubj, const int nParticles);
    ~ParticleSystem() {};

    void SetPositions(OptimizerParametersType* params);
    void GetPositions(OptimizerParametersType* params);
    void UpdateForce();
    void UpdateConstraint();
    void Integrate();

    // functor for integration
    void operator()(const VNLVector &x, VNLVector& dxdt, const double t);

    // functor for observer
    void operator()(const VNLVector &x, const double t);

private:
    const int m_nDim;
    const int m_nSubject;
    const int m_nParticles;
    const int m_nParams;
    double m_Cutoff;
    double m_Sigma2;
    double m_Viscosity;
    double m_Mu;
    int m_TimeStep;

    VNLVector m_Status;
    VNLMatrixRef m_Pos;
    VNLMatrixRef m_Vel;
    VNLMatrix m_Force;
    VNLMatrixArray m_PosArray;
};

#endif /* defined(__ParticlesGUI__myParticleDynamics__) */
