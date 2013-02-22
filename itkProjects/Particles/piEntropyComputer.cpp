//
//  piEntropyComputer.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 2/21/13.
//
//

#include "piEntropyComputer.h"
#include "iostream"

using namespace std;

namespace pi {
    void EntropyComputer::MoveToCenter() {
        int ndatasize = nsamples * nelems;
        mean.set_size(ndatasize);
        mean.fill(0);

        DataIterator3 meanIter(mean.data_block(), nsamples, nelems);
        dataIter.FirstData();
        for (int i = 0; i < ndata; i++) {
            for (int j = 0; j < nsamples; j++) {
                for (int k = 0; k < nelems; k++) {
                    meanIter.sample[k] += (dataIter.sample[k]);
                }
                meanIter.NextSample();
                dataIter.NextSample();
            }
            dataIter.NextData();
            meanIter.FirstData();
        }
        mean /= ndata;
        dataIter.FirstData();
        for (int i = 0; i < ndata; i++) {
            for (int j = 0; j < nsamples; j++) {
                for (int k = 0; k < nelems; k++) {
                    dataIter.sample[k] -= meanIter.sample[k];
                }
                meanIter.NextSample();
                dataIter.NextSample();
            }
            dataIter.NextData();
            meanIter.FirstData();
        }
    }

    // in dual space
    bool EntropyComputer::ComputeCovariance(double alpha) {
        if (ndata > nsamples * nelems) {
            cout << "not yet implemented" << endl;
            return false;
        }

        // compute DD'
        int datasize = nsamples * nelems;
        covariance.set_size(ndata, ndata);
        covariance.fill(0);

        DataIterator3 iter1(dataMatrix.data_block(), nsamples, nelems);
        DataIterator3 iter2(dataMatrix.data_block(), nsamples, nelems);

        for (int i = 0; i < ndata; i++) {
            iter2.FirstData();
            for (int j = 0; j < ndata; j++) {
                for (int k = 0; k < datasize; k++) {
                    covariance[i][j] += iter1.data[k] * iter2.data[k];
                }
                iter2.NextData();
            }
            iter1.NextData();
        }

        for (int i = 0; i < ndata; i++) {
            for (int j = 0; j < ndata; j++) {
                covariance[i][j] /= (datasize-1);
            }
            covariance[i][i] += alpha;
        }
        inverseCovariance = vnl_matrix_inverse<double>(covariance);
        return true;
    }




    bool EntropyComputer::ComputeGradient(VNLDoubleMatrix& gradient) {
        // Y*(K+I)^-1
        dataIter.FirstData();
        gradient.set_size(ndata, ndatasize);
        gradient.fill(0);
        DataIterator3 out(gradient.data_block(), ndatasize, 1);
        DataIterator3 in(dataMatrix.data_block(), ndatasize, 1);
        // start from out column

        out.FirstData();
        for (int j = 0; j < ndata; j++) {
            for (int i = 0; i < ndatasize; i++) {
                in.FirstData();
                for (int k = 0; k < ndata; k++) {
                    out.data[i] += in.data[i] * inverseCovariance[k][j];
                    in.NextData();
                }
            }
            out.NextData();
        }
        return true;
    }
}