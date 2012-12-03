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
#include "myEventCallback.h"

class ImageParticlesAlgorithm;

class ParticleSystem {
    friend class ParticleSystemObserver;

public:
    ParticleSystem(const int nSubj, const int nParticles);
    ~ParticleSystem() {};

    void SetPositions(OptimizerParametersType* params);
    void GetPositions(OptimizerParametersType* params);
    void SetConstraint(myImplicitSurfaceConstraint* constraint);
    void SetHistoryVector(VNLVectorArray* systemHistory);
    void SetCostHistoryVector(STDDoubleArray* costHistory);
    void SetEventCallback(EventCallback* callback);
    void SetContext(ImageParticlesAlgorithm* ctx);

    void Integrate();

    // functor for integration
    void operator()(const VNLVector &x, VNLVector& dxdt, const double t);

    // functor for observer
    // moved to observer
    void operator()(const VNLVector &x, const double t);

private:
    ParticleSystem() : m_nDim(0), m_nSubjects(0), m_nParticles(0), m_nParams(0) {};

    // estimate rigid transformation
    void EstimateRigidTransform(VNLMatrixRef& gPos, VNLMatrixArray& transforms, VNLMatrixArray& jacobians);

    void ApplyMatrixOperation(const double* posIn, const VNLMatrix& matrix, double* posOut);
    
    // apply global transform
    // compute forces between particles
    void UpdateSurfaceForce(VNLMatrixRef& pos, VNLMatrixRef& vel, VNLMatrix& force);

    // test code for gravity physics
    void UpdateGravityForce(VNLMatrixRef& pos, VNLMatrixRef& vel, VNLMatrix& force);

    // apply ensemble constraint
    void UpdateEnsembleForce(VNLMatrixRef& gPos, VNLMatrixRef& gVel, VNLMatrix& gForce);

    // apply constraint on implicit boundaries
    void ApplyBoundaryConditions(const VNLVector &x, const VNLMatrix& gForce, VNLMatrixRef& dpdt, VNLMatrixRef& dvdt);

    const int m_nDim;
    const int m_nSubjects;
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
    EventCallback* m_Callback;

    ImageParticlesAlgorithm* m_Context;
    myImplicitSurfaceConstraint* m_Constraint;
    VNLVectorArray* m_StatusHistory;
    STDDoubleArray* m_CostHistory;
};


#endif /* defined(__ParticlesGUI__myParticleDynamics__) */