//
//  piPlutoCore.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 5/2/13.
//
//

#include "piPlutoCore.h"

namespace pi {
    ImageIO<RealImage3> __real3IO;
    ImageIO<RealImage2> __real2IO;
    
    PlutoCore::PlutoCore(QObject* parent): QObject(parent) {

    }

    PlutoCore::~PlutoCore() {

    }

    void PlutoCore::addImage(RealImage2::Pointer images) {

    }

    void PlutoCore::initParticles(pi::Particle *particles) {

    }

    void PlutoCore::run() {
        
    }
}