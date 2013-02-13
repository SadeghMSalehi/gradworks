//
//  myImageProcessing.h
//  ParticlesGUI
//
//  Created by Joohwi Lee on 2/11/13.
//
//

#ifndef __ParticlesGUI__myImageProcessing__
#define __ParticlesGUI__myImageProcessing__

#include <iostream>
#include "piImageDef.h"

namespace pi {
    class ImageProcessing {
    public:
        // anti-aliasing, connected component, and closing morphology
        LabelImage::Pointer SmoothLabelMap(LabelImage::Pointer img);
        LabelImage::Pointer ErodedBorder(LabelImage::Pointer img);
        VectorImage::Pointer ComputeNormal(LabelImage::Pointer img);
        LabelImage::Pointer Ellipse(int* outputSize, double* center, double* radius);
        VectorImage::Pointer DistanceMap(LabelImage::Pointer img);
    };
}
#endif /* defined(__ParticlesGUI__myImageProcessing__) */
