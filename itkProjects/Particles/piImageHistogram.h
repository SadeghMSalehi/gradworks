//
//  piImageHistogram.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/20/13.
//
//

#ifndef __ParticleGuidedRegistration__piImageHistogram__
#define __ParticleGuidedRegistration__piImageHistogram__

#include <iostream>
#include "piImageDef.h"

namespace pi {
    class ImageHistogram {
    public:
        std::vector<int> binData;

    public:
        ImageHistogram(): binCount(100) {}
        
        void SetDataRange(DataReal min, DataReal max) {
            _dataMin = min;
            _dataMax = max;
        }
        void SetBinCount(int n) {
            binCount = n;
            binData.resize(binCount);
        }
        void Compute(RealImage::Pointer img);

        DataReal dataMin() { return _dataMin; }
        DataReal dataMax() { return _dataMax; }

    private:
        int binCount;
        DataReal _dataMin;
        DataReal _dataMax;
    };
}
#endif /* defined(__ParticleGuidedRegistration__piImageHistogram__) */
