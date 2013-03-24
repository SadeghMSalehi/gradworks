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
#include "itkStatisticsImageFilter.h"

namespace pi {
    template <class T>
    class ImageHistogram {
    public:
        typedef typename T::PixelType PixelType;
        std::vector<int> binData;
        PixelType dataMin, dataMax;
        PixelType rangeMin, rangeMax;

    private:
        int _binCount;
        typename T::Pointer _img;

    public:
        ImageHistogram(): _binCount(100) {
            binData.resize(_binCount);
        }

        void SetImage(typename T::Pointer img) {
            _img = img;
            typedef itk::StatisticsImageFilter<T> StatFilter;
            typename StatFilter::Pointer statFilter = StatFilter::New();
            statFilter->SetInput(_img);
            statFilter->Update();
            rangeMin = dataMin = statFilter->GetMinimum();
            rangeMax = dataMax = statFilter->GetMaximum();
        }
        void SetRange(PixelType min, PixelType max) {
            rangeMin = min;
            rangeMax = max;
        }
        void SetBinCount(int n) {
            _binCount = n;
            binData.resize(_binCount);
        }
        void Compute() {
            const PixelType* buff = _img->GetBufferPointer();
            int nSize = _img->GetPixelContainer()->Size();
            const PixelType dataMin = this->rangeMin;
            const PixelType dataRange = this->rangeMax - dataMin;
            const int binCount = this->_binCount;
            for (int i = 0; i < nSize; i++) {
                int kIdx = (buff[i]-dataMin)/(dataRange) * (binCount-1);
                binData[kIdx] ++;
            }
        }
    };
}
#endif /* defined(__ParticleGuidedRegistration__piImageHistogram__) */
