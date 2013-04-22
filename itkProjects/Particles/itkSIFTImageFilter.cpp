//
//  itkSIFTImageFilter.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/21/13.
//
//

#include "itkSIFTImageFilter.h"
#include <vnl/algo/vnl_svd.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>

using namespace std;

namespace itk {
    typedef unsigned int uint;
    typedef unsigned char uchar;
    
    void SIFTImagePCAComputer::computePCA(SIFTImage *image) {
        SIFTFeature* feature = image->GetBufferPointer();
        const uint size = image->GetPixelContainer()->Size();
        _mean.fill(0);
        const uint featureSize = _mean.size();
        for (uint i = 0; i < size; i++) {
            for (uint j = 0; j < featureSize; j++) {
                _mean[j] += feature[i][j];
            }
        }
        _mean /= size;
        
        _siftImage.set_size(featureSize, size);
        for (uint j = 0; j < featureSize; j++) {
            for (uint i = 0; i < size; i++) {
                _siftImage[j][i] = feature[i][j] - _mean[j];
            }
        }
        _siftImageCov = _siftImage * _siftImage.transpose();
        
        vnl_symmetric_eigensystem<double> siftEigen(_siftImageCov);
        cout << siftEigen.D << endl;
    }
    
}