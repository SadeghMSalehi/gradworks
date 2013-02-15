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

namespace pi {
    class ParticleSystemSolver {
        void Preprocessing(ParticleSystem& system);
        void Run(ParticleSystem& system);
    private:
        double m_t0;
        double m_dt;
        double m_t1;
    };
}
#endif /* defined(__ParticlesGUI__piParticleSystemSolver__) */
