//
//  piImageEntropy.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/12/13.
//
//

#ifndef __ParticleGuidedRegistration__piImageEntropy__
#define __ParticleGuidedRegistration__piImageEntropy__

#include <iostream>

#include "piImageDef.h"
#include "piEntropyComputer.h"

namespace pi {
    class ImageHolder {
    public:
        ImageHolder() {
            isNormalized = false;
        }
        bool isNormalized;
        RealImage::Pointer referencePatch;
        RealImage::Pointer sourceImage;
        RealImage::Pointer output;
    };

    class ImageEntropyComputer {
    public:
        ImageEntropyComputer();

        inline void operator()() {
            computeEntropy();
        }

        inline void setNormalized(bool val) {
            _isNormalized = val;
        }

        void setSize(int mData, int nSamples);
        void addSample(DataReal* sample);
        void clear();

        double entropyValue();
        void computeEntropy();

        static RealImage::Pointer computeEntropy(RealImage::Pointer, RealImage::RegionType, RealImage::Pointer, RealImage::RegionType);
        static RealImage::Pointer computeNormalizedEntropy(RealImage::Pointer, RealImage::RegionType, RealImage::Pointer, RealImage::RegionType);
        static RealImage::Pointer computeMeanSquares(RealImage::Pointer, RealImage::RegionType, RealImage::Pointer, RealImage::RegionType);
        static RealImage::Pointer computeCrossCorrelation(RealImage::Pointer, RealImage::RegionType, RealImage::Pointer, RealImage::RegionType);
    private:
        std::vector<DataReal*> _data;
        double _value;
        int _m;
        int _n;
        bool _isNormalized;

        vnl_vector<double> _meanStore;
        vnl_matrix<double> _dataStore;
        vnl_matrix<double> _covStore;
        vnl_matrix<double> _V;
        vnl_vector<double> _D;
    };
}
#endif /* defined(__ParticleGuidedRegistration__piImageEntropy__) */
