//
//  piImageRegistration.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/18/13.
//
//

#ifndef __ParticleGuidedRegistration__piImageRegistration__
#define __ParticleGuidedRegistration__piImageRegistration__

#include <iostream>
#include "piImageDef.h"

namespace pi {
    RealImage::Pointer bsplineRegistration(RealImage::Pointer srcImg, RealImage::Pointer dstImg);
}

#endif /* defined(__ParticleGuidedRegistration__piImageRegistration__) */
