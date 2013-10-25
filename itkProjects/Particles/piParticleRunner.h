//
//  piParticleRunner.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 10/25/13.
//
//

#ifndef __ParticleGuidedRegistration__piParticleRunner__
#define __ParticleGuidedRegistration__piParticleRunner__

#include <iostream>

#include "piOptions.h"
#include "piConfigFile.h"


namespace pi {
    class ParticleRunner {
    public:
        void main(Options& opts, StringVector& args);

    private:
        void print();
        void initialize(Options& opts, StringVector& args);
        void buildPatches(libconfig::Setting& setting);

        // optimization functions
        void computePatchMapping();

        ConfigFile _config;
    };
}

#endif /* defined(__ParticleGuidedRegistration__piParticleRunner__) */
