//
//  piParticleSystemSolver.h
//  ParticlesGUI
//
//  Created by Joohwi Lee on 2/15/13.
//
//

#ifndef __ParticlesGUI__piParticleSystemSolver__
#define __ParticlesGUI__piParticleSystemSolver__

#include <iostream>
#include "piParticleCore.h"
#include "piOptions.h"
#include "piParticleTrace.h"
#include "piParticleForces.h"
#include "piParticleCollision.h"
#include "piTimer.h"

namespace pi {
    class ParticleSystemSolver {
    public:
        ParticleSystemSolver();
        
        bool LoadConfig(std::istream& is);
        bool LoadConfig(const char* name);
        bool SaveConfig(const char* name);

        // load parameters for execution parameters not for images or labels
        bool LoadParameters(std::istream& is);

        void Preprocessing();

        void PreprocessingWithLabels();

        // convenient function with only internal forces
        void SpreadParticles();
        void SpreadMultiParticles();

        void SetupCollisionHandlers();
        
        void Run();
        void Setup();
        void SetupParameters();

        void RunLoopBegin();
        int  RunStep();
        void RunLoopEnd();
        void Stop();

        void PrintPoints();
        
        pi::RealImage::Pointer WarpImage(int i, int j = -1);
        pi::LabelImage::Pointer WarpLabel(int i, int j = -1);


    private:
        void RunStepBegin();
        void RunStepEnd();
        void CheckEnergy();

    public:
        Options m_Options;
        ParticleSystem m_System;

        EntropyInternalForce internalForce;
        EnsembleForce ensembleForce;
        IntensityForce intensityForce;

        ParticleTrace trace;
        std::vector<ParticleCollision> collisionHandlers;
        std::vector<ParticleMultiCollision> multiCollisionHandlers;
        Timer timer;

        // system state store
        std::string systemSnapshot;
        std::string traceFile;

        // options
        bool useEnsemble;
        bool useIntensity;
        bool noInternal;
        bool noBoundary;
        bool traceOn;
        bool verbose;
        bool useVelocity;
        bool useAlignment;
        bool continueToRun;

        // time management
        DataReal t0;
        DataReal t1;
        DataReal dt;
        DataReal t;

    private:
        DataReal _previousImageEnergy;
    };
}
#endif /* defined(__ParticlesGUI__piParticleSystemSolver__) */
