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
#include "itkRescaleIntensityImageFilter.h"

class vtkPolyData;

namespace pi {
    class ImageProcessing {
    public:
        // anti-aliasing, connected component, and closing morphology
        LabelImage::Pointer SmoothLabelMap(LabelImage::Pointer img);
        LabelImage::Pointer ErodedBorder(LabelImage::Pointer img);
        
        GradientImage::Pointer ComputeGaussianGradient(LabelImage::Pointer img, double sigma = -1);
        GradientImage::Pointer ComputeGradient(LabelImage::Pointer img);
        GradientImage::Pointer ComputeGaussianGradient(DoubleImage::Pointer img, double sigma = -1);
        GradientImage::Pointer ComputeGradient(DoubleImage::Pointer img);
        DoubleImage::Pointer ComputeMagnitudeMap(VectorImage::Pointer img);
        DoubleImage::Pointer ComputeMagnitudeMap(GradientImage::Pointer img);

        DoubleImage::Pointer ComputeGaussianGradientMagnitude(DoubleImage::Pointer img, double sigma = -1);
        LabelImage::Pointer Ellipse(int* outputSize, double* center, double* radius);
        VectorImage::Pointer DistanceMap(LabelImage::Pointer img);
        vtkPolyData* ConvertToMesh(LabelImage::Pointer image);
        DoubleImage::Pointer NormalizeIntensity(DoubleImage::Pointer image, LabelImage::Pointer label);

        template <class T>
        typename T::Pointer RescaleIntensity(typename T::Pointer srcImg, typename T::PixelType min, typename T::PixelType max) const {
            typedef itk::RescaleIntensityImageFilter<T> RescaleFilter;
            typename RescaleFilter::Pointer filter = RescaleFilter::New();
            filter->SetInput(srcImg);
            filter->SetOutputMinimum(min);
            filter->SetOutputMaximum(max);
            filter->Update();
            return filter->GetOutput();
        }
    };
}
#endif /* defined(__ParticlesGUI__myImageProcessing__) */
