//
//  piParticleTools.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 2/16/13.
//
//

#ifndef __ParticleGuidedRegistration__piParticleTools__
#define __ParticleGuidedRegistration__piParticleTools__

#include <iostream>
#include "piImageDef.h"
#include "piParticleCore.h"

namespace pi {
    template <class T>
    void MarkAtImage(T& data, int n, LabelImage::Pointer image, LabelPixel val) {
        for (int i = 0; i < n; i++) {
            Particle& pi = data[i];
            IntIndex idx;
            fordim (k) {
                idx[k] = pi.x[k] + 0.5;
            }
            image->SetPixel(idx, val);
        }
    }
}
#endif /* defined(__ParticleGuidedRegistration__piParticleTools__) */
