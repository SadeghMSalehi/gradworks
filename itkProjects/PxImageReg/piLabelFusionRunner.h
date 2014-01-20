//
//  piLabelFusionRunner.h
//  PxImageReg
//
//  Created by Joohwi Lee on 1/16/14.
//
//

#ifndef __PxImageReg__piLabelFusionRunner__
#define __PxImageReg__piLabelFusionRunner__

#include <iostream>
#include "piImageDef.h"

namespace pi {
    /// @brief Perform majority voting for the final segmentation
    /// @param labels A vector of label images
    /// @param images A vector of intensity images
    /// @param targetImage the target intensity image to be segmented
    LabelImage::Pointer performLabelFusion(LabelImageVector& labels, RealImageVector& images, RealImage::Pointer targetImage);
};

#endif /* defined(__PxImageReg__piLabelFusionRunner__) */
