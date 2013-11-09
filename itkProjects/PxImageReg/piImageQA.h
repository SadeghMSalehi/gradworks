//
//  piImageQA.h
//  PxImageReg
//
//  Created by Joohwi Lee on 11/9/13.
//
//

#ifndef __PxImageReg__piImageQA__
#define __PxImageReg__piImageQA__

#include <iostream>
#include <vector>

#include "piImageIO.h"
#include "piImageDef.h"

namespace pi {
    class ImageQA {
    public:
        int axis;
        int slice;
        float alpha;

        /// load settings for alpha value, label colors, etc
        void loadConfig(std::string file);
        void load(std::string image, std::string label);
        void save(std::string output);

    private:
        RealImage3::Pointer image;
        LabelImage3::Pointer label;
    };
}
#endif /* defined(__PxImageReg__piImageQA__) */
