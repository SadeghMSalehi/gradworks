//
//  itkSIFTImageFilter.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/21/13.
//
//

#ifndef __ParticleGuidedRegistration__itkSIFTImageFilter__
#define __ParticleGuidedRegistration__itkSIFTImageFilter__

#include <iostream>
#include <itkImageToImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkVectorLinearInterpolateImageFunction.h>
#include "piImageDef.h"
#include "piImageEntropyComputer.h"
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>

using namespace pi;

namespace itk {

    typedef itk::Vector<u_int8_t,128> SIFTFeature;
    typedef itk::Image<SIFTFeature,DIMENSIONS> SIFTImage;

    class SIFTImageFilter: public ImageToImageFilter<GradientImage,SIFTImage> {
    public:
        typedef SIFTImageFilter Self;
        typedef ImageToImageFilter<GradientImage,SIFTImage> Superclass;
        typedef SmartPointer<Self> Pointer;
        typedef SmartPointer<const Self> ConstPointer;

        itkNewMacro(SIFTImageFilter);
        itkTypeMacro(SIFTImageFilter,ImageToImageFilter);
        
    protected:
        typedef GradientImage::RegionType RegionType;
        typedef GradientImage::IndexType IndexType;

        SIFTImageFilter() {

        }

        virtual ~SIFTImageFilter() {

        }

        /** Does the real work. */
        virtual void ThreadedGenerateData(const typename Superclass::OutputImageRegionType& outputRegionForThread, ThreadIdType threadId) {
            GradientImage::ConstPointer input = this->GetInput();
            SIFTImage::Pointer output = this->GetOutput();

            GradientInterpolatorType::Pointer pixelSampler = GradientInterpolatorType::New();
            pixelSampler->SetInputImage(input);

            ImageRegionIteratorWithIndex<SIFTImage> outIter(output, outputRegionForThread);
            RegionType siftWindow;
            siftWindow.SetSize(0, 16);
            siftWindow.SetSize(1, 16);

            ImageGradientHistogram siftHisto;
            for (outIter.GoToBegin(); !outIter.IsAtEnd(); ++outIter) {
                IndexType idx = outIter.GetIndex();
                siftWindow.SetIndex(0, (idx[0]-8<0) ? 0 : idx[0]-8);
                siftWindow.SetIndex(1, (idx[1]-8<0) ? 0 : idx[1]-8);

                SIFTImage::PixelType feature;
                feature.Fill(0);

                siftHisto.setRegion(siftWindow);
                siftHisto.computeHistogram(input);

                memcpy(feature.GetDataPointer(), siftHisto.histogram(), 128);
                outIter.Set(feature);
            }
        }

    private:
        SIFTImageFilter(const SIFTImageFilter&);
        void operator=(const SIFTImageFilter&);
    };
    
    class SIFTImagePCAComputer {
    public:
        void computePCA(SIFTImage* image);
        
    private:
        vnl_vector<double> _mean;
        vnl_matrix<double> _siftImage;
        vnl_matrix<double> _siftImageCov;
    };
}
#endif /* defined(__ParticleGuidedRegistration__itkSIFTImageFilter__) */
