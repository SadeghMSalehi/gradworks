//
//  ParticleWarp.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 10/31/13.
//
//

#ifndef __ParticleGuidedRegistration__ParticleWarp__
#define __ParticleGuidedRegistration__ParticleWarp__

#include <iostream>
#include "piConfigFile.h"
#include "piParticleRunner.h"

namespace pi {
    class ParticleMesh {
    public:
        void constructNeighbors(int regionId, int nPx, PxSubj& subj, double cutoff, PxGlobal::Neighbors& neighbors);
    };

    class ParticleWarp {
    public:
        int controlSpacing;
        
        LabelImage::Pointer reference;
        DisplacementFieldType::Pointer displacementField;

        ParticleWarp(): controlSpacing(4) {}

        void setParameters(ConfigFile& config);
        void estimateBsplineWarp(Px::Vector& src, Px::Vector& dst);
        LabelImage::Pointer warpLabel(LabelImage::Pointer input);
        RealImage::Pointer warpImage(RealImage::Pointer input);

    private:
        DisplacementFieldPointSetType::Pointer m_FieldPoints;
        FieldTransformType::Pointer transform;
    };
}
#endif /* defined(__ParticleGuidedRegistration__ParticleWarp__) */
