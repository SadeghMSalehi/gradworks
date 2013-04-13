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

        void setSize(int mData, int nSamples);
        void addSample(DataReal* sample);
        void clear();
        
        double entropyValue();
        void computeEntropy();

        static RealImage::Pointer computeEntropy(RealImage::Pointer, RealImage::RegionType, RealImage::Pointer, RealImage::RegionType);

    private:
        std::vector<DataReal*> _data;
        double _value;
        int _m;
        int _n;

        vnl_vector<DataReal> _meanStore;
        vnl_matrix<DataReal> _dataStore;
        vnl_matrix<DataReal> _covStore;
        vnl_matrix<DataReal> _V;
        vnl_vector<DataReal> _D;
    };
}
#endif /* defined(__ParticleGuidedRegistration__piImageEntropy__) */
