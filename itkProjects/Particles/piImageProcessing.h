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

class vtkPolyData;

namespace pi {
    class ImageProcessing {
    public:
        // anti-aliasing, connected component, and closing morphology
        LabelImage::Pointer SmoothLabelMap(LabelImage::Pointer img);
        LabelImage::Pointer ErodedBorder(LabelImage::Pointer img);
        VectorImage::Pointer ComputeGradient(LabelImage::Pointer img);
        LabelImage::Pointer Ellipse(int* outputSize, double* center, double* radius);
        VectorImage::Pointer DistanceMap(LabelImage::Pointer img);
        DoubleImage::Pointer ComputeMagnitudeMap(VectorImage::Pointer img);
        vtkPolyData* ConvertToMesh(LabelImage::Pointer image);
        DoubleImage::Pointer NormalizeIntensity(DoubleImage::Pointer image, LabelImage::Pointer label);
    };
}
#endif /* defined(__ParticlesGUI__myImageProcessing__) */
