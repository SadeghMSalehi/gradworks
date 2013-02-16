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

namespace pi {
    class ParticleSystemSolver {
    public:
        bool LoadConfig(const char* name);
        bool SaveConfig(const char* name);
        void Preprocessing();
        void Run();
        
        ImageContext& GetImageContext() {
            return m_ImageContext;
        }

        Options m_Options;
        ParticleSystem m_System;
        ImageContext m_ImageContext;
    };
}
#endif /* defined(__ParticlesGUI__piParticleSystemSolver__) */
