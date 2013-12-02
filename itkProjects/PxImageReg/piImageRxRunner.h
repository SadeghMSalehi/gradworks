//
//  piImageRxRunner.h
//  PxImageReg
//
//  Created by Joohwi Lee on 11/30/13.
//
//

#ifndef __PxImageReg__piImageRxRunner__
#define __PxImageReg__piImageRxRunner__

#include <iostream>
#include "piOptions.h"


namespace pi {
    class ImageRx {
    public:
        ImageRx(Options& opts, StringVector& args) {}

        // main entry functions
        void mainRigidRegistration(Options& opts, StringVector& args);
    };

}
#endif /* defined(__PxImageReg__piImageRxRunner__) */
