//
//  piParticleTrace.h
//  ParticlesGUI
//
//  Created by Joohwi Lee on 2/16/13.
//
//

#ifndef __ParticlesGUI__piParticleTrace__
#define __ParticlesGUI__piParticleTrace__

#include <iostream>
#include "piParticleCore.h"

namespace pi {
    class ParticleTrace {
    public:
        void Add(double t, ParticleArray& array);
        void Write(std::ostream& os);
        void Read(std::istream& is);
        void Read(std::istream& is, ParticleVector& trace);
    private:
        ParticleVector m_Trace;
    };
}

#endif /* defined(__ParticlesGUI__piParticleTrace__) */
