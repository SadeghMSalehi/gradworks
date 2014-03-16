//
//  piDemonsRunner.h
//  PxImageReg
//
//  Created by Joohwi Lee on 11/20/13.
//
//

#ifndef __PxImageReg__piDemonsRunner__
#define __PxImageReg__piDemonsRunner__

#include <iostream>
#include "piOptions.h"
#include "piImageDef.h"
#include "libconfig.h++"

namespace pi {
    class DemonsRunner {
    public:
        DemonsRunner(Options& opts, StringVector& args);

        void buildPatches(libconfig::Setting& setting);
        void computePatchMapping();
        void computeOpticalFlowMapping();

        DisplacementFieldType::Pointer computeDemonsFlow(RealImage::Pointer fixed, RealImage::Pointer moving, double dt);

        RealImage::Pointer deformImage(RealImage::Pointer input, DisplacementFieldType::Pointer displacement, RealImage::Pointer refImage);

        /// Warp a vector image in according to the displacement image
        /// The displacement image defines a mapping from the output image coordinates to the input image coordinates.
        /// The mapping is linear so that x_i = x_o + d(x_o)
        /// The output pixel is resampled from the corresponding input point: p(x_o) = p(x_o + d(x_o)).
        DisplacementFieldType::Pointer deformVectorImage(
            DisplacementFieldType::Pointer input,
            DisplacementFieldType::Pointer displacement, RealImage::Pointer refImage);

        /// Compute an intermediate image to record the transformed coordinate
        ///
        /// \param sourceImage an image defines the space of the fixed image
        ///
        DisplacementFieldType::Pointer computeCoordImage(RealImage::Pointer sourceImage);

        /// Deform an image with transformed coordinates
        ///
        /// \param sourceImage an image to be warped
        /// \param coordImage an image that contains point coordinates to sample
        RealImage::Pointer resampleImage(RealImage::Pointer sourceImage, DisplacementFieldType::Pointer coordImage);

        DisplacementFieldType::Pointer resampleField(DisplacementFieldType::Pointer currentField, DisplacementFieldType::Pointer resamplingField);

    private:
        libconfig::Config _config;
    };
}

#endif /* defined(__PxImageReg__piDemonsRunner__) */
