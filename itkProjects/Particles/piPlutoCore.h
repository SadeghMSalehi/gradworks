//
//  piPlutoCore.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 5/2/13.
//
//

#ifndef __ParticleGuidedRegistration__piPlutoCore__
#define __ParticleGuidedRegistration__piPlutoCore__

#include <iostream>
#include <vector>
#include "piParticle.h"
#include "piImageDef.h"

namespace pi {
    class PlutoCore: public QObject {
        Q_OBJECT
        
    public:
        PlutoCore(QObject* parent = NULL);
        virtual ~PlutoCore();
        
        void addImage(RealImage::Pointer images);
        void initParticles(Particles*);
        void run();
        
    private:
        std::vector<RealImage::Pointer> _images;

    };
}
#endif /* defined(__ParticleGuidedRegistration__piPlutoCore__) */
