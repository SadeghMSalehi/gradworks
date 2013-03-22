//
//  piImageHistogram.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/20/13.
//
//

#include "piImageHistogram.h"

namespace pi {
    void ImageHistogram::Compute(RealImage::Pointer img) {
        const DataReal* buff = img->GetBufferPointer();
        int nSize = img->GetPixelContainer()->Size();
        const DataReal dataMin = this->_dataMin;
        const DataReal dataRange = _dataMax - dataMin;
        const int binCount = this->binCount;
        binData.resize(binCount);
        for (int i = 0; i < nSize; i++) {
            int kIdx = (buff[i]-dataMin)/(dataRange) * binCount;
            binData[kIdx] ++;
        }
    }
}