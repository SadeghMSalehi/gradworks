//
//  piImageAlign.h
//  PxImageReg
//
//  Created by Joohwi Lee on 11/17/13.
//
//

#ifndef __PxImageReg__piImageAlign__
#define __PxImageReg__piImageAlign__

#include <iostream>
#include "piOptions.h"

/***
This class is written to provide a gradient-based affine alignment for multi-modal images. First, this extracts a gradient image then perform the gradient's CDF based alignment.
 ***/

namespace pi {
    class ImageAlign {
        ImageAlign();

        void computeGradientMagnitude(Options& opts, StringVector& args);
        void setupOptions(Options& opts);
        void main(Options& opts, StringVector& args);
    };
}

#endif /* defined(__PxImageReg__piImageAlign__) */
