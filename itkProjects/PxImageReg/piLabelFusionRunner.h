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
    LabelImage::Pointer performLabelFusion(LabelImageVector& labels, RealImageVector& images, RealImage::Pointer targetImage);
};

#endif /* defined(__PxImageReg__piLabelFusionRunner__) */
