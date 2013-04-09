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
        void Preprocessing();

        void Run();
        void Setup();
        void RunStepBegin();
        int  RunStep();
        void RunStepEnd();

        void PrintPoints();
        
        pi::RealImage::Pointer WarpImage(int i, int j = -1);
        pi::LabelImage::Pointer WarpLabel(int i, int j = -1);

        
        ImageContext& GetImageContext() {
            return m_ImageContext;
        }

        Options m_Options;
        ParticleSystem m_System;
        ImageContext m_ImageContext;

        EntropyInternalForce internalForce;
        EnsembleForce ensembleForce;
        IntensityForce intensityForce;

        ParticleTrace trace;
        std::vector<ParticleCollision> collisionHandlers;
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

        // time management
        DataReal t0;
        DataReal t1;
        DataReal dt;
        DataReal t;
    };
}
#endif /* defined(__ParticlesGUI__piParticleSystemSolver__) */
