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
        PixelType pixelMean, pixelVariance, pixelStdev;
        double fitPercentile;

    private:
        int _binCount;
        typename T::Pointer _img;

    public:
        ImageHistogram(): _binCount(100) {
            fitPercentile = 0.01;
            binData.resize(_binCount);
        }
        void FitRange() {
            // fit rangeMin and rangeMax to 1% from dataMin and dataMax
            const int binCount = 10000;
            std::vector<int> binData(binCount);
            PixelType* buff = _img->GetBufferPointer();
            const int nSize = _img->GetPixelContainer()->Size();
            const double dataMin = this->dataMin;
            const double dataRange = (this->dataMax - this->dataMin);
            for (int i = 0; i < nSize; i++) {
                unsigned short bin = (unsigned short)(binCount*(double(*buff) - dataMin)/(dataRange));
                binData[bin]++;
                buff++;
            }

            int marginCount = nSize * fitPercentile;
            int acc = 0;
            for (int i = 0; i < binCount; i++) {
                acc += binData[i];
                if (acc > marginCount) {
                    rangeMin = double(i)/double(binCount)*(dataRange)+dataMin;
                    break;
                }
            }
            acc = 0;
            for (int i = binCount - 1; i >= 0; i--) {
                acc += binData[i];
                if (acc > marginCount) {
                    rangeMax = double(i+1)/double(binCount)*(dataRange)+dataMin;
                    break;
                }
            }
        }
        void SetImage(typename T::Pointer img, bool fitRange = true) {
            _img = img;
            if (fitRange) {
                typedef itk::StatisticsImageFilter<T> StatFilter;
                typename StatFilter::Pointer statFilter = StatFilter::New();
                statFilter->SetInput(_img);
                statFilter->Update();
                dataMin = statFilter->GetMinimum();
                dataMax = statFilter->GetMaximum();
                pixelMean = statFilter->GetMean();
                pixelVariance = statFilter->GetVariance();
                pixelStdev = statFilter->GetSigma();
                FitRange();
            }
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
        inline PixelType NormalizePixel(PixelType p) {
            return (p-pixelMean)/(pixelStdev);
        }
    };
}
#endif /* defined(__ParticleGuidedRegistration__piImageHistogram__) */
