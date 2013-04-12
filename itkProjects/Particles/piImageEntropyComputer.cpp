//
//  piImageEntropy.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/12/13.
//
//

#include "piImageEntropyComputer.h"

namespace pi {
    ImageEntropyComputer::ImageEntropyComputer() {
        _value = 0;
        _m = _n = 0;
    }

    void ImageEntropyComputer::setSize(int mData, int nSamples) {
        _m = mData;
        _n = nSamples;

        _meanStore.set_size(_n);
        _dataStore.set_size(_m, _n);
        _covStore.set_size(_m, _m);
    }

    void ImageEntropyComputer::addSample(DataReal* sampleBuffer) {
        _data.push_back(sampleBuffer);
    }

    void ImageEntropyComputer::clear() {
        _data.clear();
        _value = 0;
    }

    double ImageEntropyComputer::entropyValue() {
        return _value;
    }

    void ImageEntropyComputer::computeEntropy() {
        const int n = _n;
        const int m = _m;

        _meanStore.fill(0);
        for (int j = 0; j < m; j++) {
            DataReal* mean = _meanStore.data_block();
            DataReal* pixels = _data[j];
            for (int i = 0; i < n; i++, mean++) {
                *mean += *pixels;
            }
        }
        _meanStore /= m;


        DataReal* data = _dataStore.data_block();
        for (int j = 0; j < m; j++) {
            DataReal* mean = _meanStore.data_block();
            DataReal* pixels = _data[j];
            for (int i = 0; i < n; i++, data++, mean++) {
                *data = *pixels - *mean;
            }
        }

        DataReal* cov = _covStore.data_block();
        // i-th row
        for (int i = 0; i < _m; i++) {
            // k-th column
            for (int k = 0; k < _m; k++) {
                // j-th element
                for (int j = 0; j < _n; j++) {
                    *cov = _dataStore[i][j] * _dataStore[k][j];
                }
                ++cov;
            }
        }

        vnl_symmetric_eigensystem_compute(_covStore, _V, _D);
        _value = _D.sum();
    }
}