//
//  piImageAlign.cpp
//  PxImageReg
//
//  Created by Joohwi Lee on 11/17/13.
//
//

#include "piImageAlign.h"

namespace pi {
    ImageAlign::ImageAlign() {

    }

    void ImageAlign::setupOptions(pi::Options &opts) {
        opts.addOption("--gradmag", "compute a gradient magnitude image", SO_NONE);
        opts.addOption("--sigma", "sigma to perform a gaussian smoothing", SO_REQ_SEP);
    }

    void ImageAlign::main(Options& opts, StringVector& args) {
        if (opts.GetBool("--gradmag")) {
            computeGradientMagnitude(opts, args);
        }
    }
}