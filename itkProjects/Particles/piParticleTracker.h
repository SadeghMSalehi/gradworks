//
//  piParticleTracker.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/30/13.
//
//

#ifndef __ParticleGuidedRegistration__piParticleTracker__
#define __ParticleGuidedRegistration__piParticleTracker__

#include <iostream>
#include <vector>
#include "piImageDef.h"
#include "piParticle.h"
#include "piImagePatch.h"

namespace pi {
    class ParticleTracker2 {
    public:
        ParticleTracker2();
        ~ParticleTracker2();

        void addImage(RealImage::Pointer image);
        void setParticles(Particle* particles);

    private:
        std::vector<RealImage::Pointer> _images;
//        ImageSamples<RealImage> _samples;
    };
}
#endif /* defined(__ParticleGuidedRegistration__piParticleTracker__) */
