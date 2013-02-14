//
//  piParticleHydroDynamics.h
//  ParticlesGUI
//
//  Created by Joohwi Lee on 2/14/13.
//
//

#ifndef __ParticlesGUI__piParticleHydroDynamics__
#define __ParticlesGUI__piParticleHydroDynamics__

#include <iostream>
#include "piParticleCore.h"

namespace pi {

    class HydroDynamics {
    public:
        void ComputeValues(ParticleSubject& subj);
    };
}
#endif /* defined(__ParticlesGUI__piParticleHydroDynamics__) */
