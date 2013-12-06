//
//  piImageProc.h
//  PxImageReg
//
//  Created by Joohwi Lee on 11/5/13.
//
//

#ifndef __PxImageReg__piImageProc__
#define __PxImageReg__piImageProc__

#include <iostream>
#include "piImageDef.h"

#include <itkSmoothingRecursiveGaussianImageFilter.h>

namespace pi {
    /// compute distance map that contains an offset to the closest point
    /// @bmap LabelImage a binary mask
    VectorImage::Pointer ComputeDistanceMap(LabelImage::Pointer bmap);

    /// compute gaussian gradient with sigma from a label image
    GradientImage::Pointer ComputeGaussianGradient(LabelImage::Pointer bmap, double sigma);

    /// compute gaussian gradient with sigma from a label image
    GradientImage::Pointer ComputeGaussianGradient(RealImage::Pointer img, double sigma);


    /// allocate three dimensional image to store multiple 2d images
    LabelImage3::Pointer CreateImage3(LabelImage::Pointer refImage, int m);

    
    /// extract a slice of a volume
    /// @param volume RealImage3
    /// @param dim dimension to slice
    RealImage2Vector SliceVolume(RealImage3::Pointer volume, int dim);

    template <typename T>
    typename T::Pointer applyGaussian(typename T::Pointer input, double sigma) {
        typedef itk::SmoothingRecursiveGaussianImageFilter<T> GaussianFilter;
        typename GaussianFilter::Pointer filter = GaussianFilter::New();
        filter->SetSigma(sigma);
        filter->SetInput(input);
        filter->Update();
        typename T::Pointer output = filter->GetOutput();
        output->DisconnectPipeline();
        return output;
    }

}
#endif /* defined(__PxImageReg__piImageProc__) */
