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
#include <QObject>

#include "piParticle.h"
#include "piImageDef.h"
#include "piImageIO.h"
#include "piImagePatch.h"

namespace pi {
    extern ImageIO<RealImage2> __real2IO;
    extern ImageIO<RealImage3> __real3IO;

    typedef ImageSamples<RealImage2, GradientImage2> RealSamples2;
    
    class PlutoCore: public QObject {
        Q_OBJECT
        
    public:
        PlutoCore(QObject* parent = NULL);
        virtual ~PlutoCore();
        
        void setImages(RealImage2Vector images);
        void setInitialParticles(ParticleVector initialParticles);
        
        void track(int model, int test);
        
        void run();
        void initialize();

        std::vector<ParticleVector>& getParticles();
        
    private:
        RealImage2Vector _images;
        GradientImage2Vector _gradientImages;
        ParticleVector _initialParticles;
        std::vector<ParticleVector> _particles;
        RealSamples2* _samples;
    };
}
#endif /* defined(__ParticleGuidedRegistration__piPlutoCore__) */
