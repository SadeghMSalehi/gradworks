//
//  piEntropyComputer.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 2/21/13.
//
//

#ifndef __ParticleGuidedRegistration__piEntropyComputer__
#define __ParticleGuidedRegistration__piEntropyComputer__

#include <iostream>
#include "piImageDef.h"


namespace pi {
    class DataIterator3 {
    public:
        double* const databegin;
        double* data;
        double* sample;
        int samplestride;
        int elemstride;

        DataIterator3(): databegin(NULL) {
            data = sample = NULL;
            samplestride = 0;
            elemstride = 0;
        }
        DataIterator3(double* begin, int ns, int ne)
        : databegin(begin), data(begin), sample(begin), samplestride(ns*ne), elemstride(ne) {}

        double& At(int i, int j, int k) {
            data = &databegin[i*samplestride];
            sample = &data[j*elemstride];
            return sample[k];
        }
        double& At(int i, int j) {
            data = &databegin[i*samplestride];
            sample = &data[j*elemstride];
            return *sample;
        }
        double* NextSample() {
            sample += elemstride;
            return sample;
        }
        double* PrevSample() {
            sample -= elemstride;
            return sample;
        }
        double* FirstSample() {
            sample = data;
            return sample;
        }
        double* PrevData() {
            data -= samplestride;
            sample = data;
            return data;
        }
        double* NextData() {
            data += samplestride;
            sample = data;
            return data;
        }
        double* FirstData() {
            data = databegin;
            sample = data;
            return data;
        }
        double* NextDataSample() {
            sample += samplestride;
            data += samplestride;
            return sample;
        }
        double* PrevDataSample() {
            sample -= samplestride;
            data -= samplestride;
            return sample;
        }
    };

    class EntropyComputer {
    public:
        VNLDoubleMatrix dataMatrix;
        DataIterator3 dataIter;

        int ndata;
        int nsamples;
        int nelems;
        int ndatasize;
        VNLDoubleVector mean;
        VNLDoubleMatrix covariance;
        VNLDoubleMatrix inverseCovariance;
        VNLDoubleMatrix gradient;

        EntropyComputer(int data, int samples, int elems)
            : dataMatrix(data,samples*elems), dataIter(dataMatrix.data_block(),samples, elems),
                ndata(data), nsamples(samples), nelems(elems), ndatasize(nsamples*nelems) {}
        ~EntropyComputer() {}

        void MoveToCenter();
        // in dual space
        bool ComputeCovariance(double alpha = 1);
        bool ComputeGradient();

    private:
        EntropyComputer() :dataIter(NULL,0,0){};
        void operator=(const EntropyComputer& ec) {};
    };
}
#endif /* defined(__ParticleGuidedRegistration__piEntropyComputer__) */
