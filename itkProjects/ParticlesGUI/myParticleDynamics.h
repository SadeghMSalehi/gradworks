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
#include "QElapsedTimer"

class ImageParticlesAlgorithm;

namespace my {
    // experiment options
    struct ParticleSystemOptions {
        ParticleSystemOptions() {
            this->applySurfaceEntropyToFirstOnly = false;
            this->applyBoundaryConditionToFirstOnly = false;
            this->useEnsembleForce = false;
            this->useParticlePhysics = true;
            this->useBoundaryCondition = true;
            this->useSurfaceForce = true;
            this->useImageForce = false;
            this->useAdaptiveSampling = false;
        }
        bool applySurfaceEntropyToFirstOnly;
        bool applyBoundaryConditionToFirstOnly;
        bool useEnsembleForce;
        bool useParticlePhysics;
        bool useBoundaryCondition;
        bool useSurfaceForce;
        bool useImageForce;
        bool useAdaptiveSampling;
    };
    
    class ParticleSystem {
//        friend class ParticleSystemObserver;
        
    public:
        ParticleSystem(const int nSubj, const int nParticles);
        ~ParticleSystem() {};
        
        void SetPositions(OptimizerParametersType* params);
        void GetPositions(OptimizerParametersType* params);
        void SetHistoryVector(VNLVectorArray* systemHistory);
        void SetCostHistoryVector(STDDoubleArray* costHistory);
        void SetEventCallback(EventCallback* callback);
        void SetContext(ImageParticlesAlgorithm* ctx);
        mySetMacro(Options, ParticleSystemOptions);
        
        
        // functor for integration
        void operator()(const VNLVector &x, VNLVector& dxdt, const double t);
        
        // functor for observer
        // moved to observer
        void operator()(const VNLVector &x, const double t);
        
        void Integrate();
        
        void Integrate(VNLVector& status, double dt, double t0, double t1, int odeMethod = 1);

    private:
        ParticleSystem() : m_nDim(0), m_nSubjects(0), m_nParticles(0), m_nParams(0),
        m_Pos(0,0,NULL),
        m_Vel(0,0,NULL),
        m_dpdt(0,0,NULL),
        m_dvdt(0,0,NULL) {};
        
        // estimate rigid transformation
        void EstimateRigidTransform(VNLMatrixRef& gPos, VNLMatrixArray& transforms, VNLMatrixArray& jacobians);
        
        void ApplyMatrixOperation(const double* posIn, const VNLMatrix& matrix, double* posOut);
        
        void UpdateTransform(VNLMatrix& tPos);
        
        // TPS, EBS transform
        void UpdateKernelTransform();
        
        // BSpline-based displacement transform
        void UpdateBSplineEnsemble();
        
        // apply global transform
        // compute forces between particles
        void UpdateSurfaceForce();
        
        // test code for gravity physics
        void UpdateGravityForce();
        
        // apply ensemble constraint
        void UpdateEnsembleForce();
        
        // apply image term
        void UpdateImageForce();

        void UpdateDraggingForce();
        
        // apply constraint on implicit boundaries
        void ApplyBoundaryConditions();


        const int m_nDim;
        const int m_nSubjects;
        const int m_nParticles;
        const int m_nParams;

        int m_ForceType;
        
        double m_Cutoff;
        double m_Sigma;
        double m_Viscosity;
        double m_Mu;
        double m_COR;
        double m_GradientScale;
        
        int m_TimeStep;
        
        
        VNLMatrixRef m_Pos;
        VNLMatrixRef m_Vel;
        VNLMatrixRef m_dpdt;
        VNLMatrixRef m_dvdt;
        
        VNLVector m_Status;
        VNLMatrix m_Force;
        EventCallback* m_Callback;
        
        ImageParticlesAlgorithm* m_Context;
        ImplicitSurfaceConstraint* m_Constraint;
        VNLVectorArray* m_StatusHistory;
        STDDoubleArray* m_CostHistory;
        
        ParticleSystemOptions m_Options;

        mutable QElapsedTimer m_Timer;
    };
}

#endif /* defined(__ParticlesGUI__myParticleDynamics__) */