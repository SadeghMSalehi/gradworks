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

#include "piImageIO.h"

using namespace std;
using namespace pi;

namespace itk {
    typedef unsigned int uint;
    typedef unsigned char uchar;
    
    void SIFTImagePCAComputer::computePCA(SIFTImage *image) {
        SIFTFeature* feature = image->GetBufferPointer();
        const uint size = image->GetPixelContainer()->Size();
        const uint featureSize = feature[0].GetNumberOfComponents();
        _mean.set_size(featureSize);
        _mean.fill(0);

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
        this->D = siftEigen.D;
        this->V = siftEigen.V;
    }

    RealImage::Pointer SIFTImagePCAComputer::computePCImage(SIFTImage* image, int k) {
        ImageIO<RealImage> io;
        RealImage::Pointer outputImg = io.NewImageT(image->GetBufferedRegion().GetSize());

        SIFTFeature* feature = image->GetBufferPointer();
        DataReal* output = outputImg->GetBufferPointer();
        const uint size = image->GetPixelContainer()->Size();
        const uint featureSize = feature[0].GetNumberOfComponents();
        for (uint i = 0; i < size; i++) {
            uchar* data = feature->GetDataPointer();
            float sum = 0;
            for (uint j = 0; j < featureSize; j++) {
                sum += this->V[j][featureSize - 1 - k] * data[j];
            }
            *output = std::abs(sum);

            output++;
            feature++;
        }
        return outputImg;
    }
}