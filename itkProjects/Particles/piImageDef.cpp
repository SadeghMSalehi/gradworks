//
//  piImageDef.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 5/2/13.
//
//

#include <itkExtractImageFilter.h>
#include "piImageDef.h"
#include "piImageIO.h"
#include "piImageHistogram.h"

namespace pi {
    typedef itk::ExtractImageFilter<RealImage3, RealImage2> FilterType;

    // external variable for easy access
    RealImageTools __realImageTools;

    RealImage3::Pointer RealImageTools::normalizeIntensity(RealImage3::Pointer image, double percentile) {
        ImageIO<RealImage3> io;
        RealImage3::Pointer output = io.CopyImage(image);

        ImageHistogram<RealImage3> histo;
        histo.fitPercentile = percentile;
        histo.SetImage(image);

        DataReal* inBuf = image->GetBufferPointer();
        DataReal* outBuf = output->GetBufferPointer();
        const int nPixels = image->GetPixelContainer()->Size();

        for (int i = 0; i < nPixels; i++) {
            outBuf[i] = (histo.rangeMax - histo.rangeMin) * inBuf[i] / (histo.dataMax - histo.dataMin) + histo.rangeMin;
        }

        ImageHistogram<RealImage3> histo2;
        histo2.SetImage(output);
        for (int i = 0; i < nPixels; i++) {
            outBuf[i] = histo2.NormalizePixel(outBuf[i]);
        }
        return output;
    }
    
    // split a volume image into slices
    RealImage2Vector RealImageTools::sliceVolume(RealImage3::Pointer volume, int dim) {
        RealImage2Vector slices;
        RealImage3::RegionType region = volume->GetBufferedRegion();
        RealImage3::SizeType sz = region.GetSize();
        const RealImage3::IndexType idx = region.GetIndex();
        
        const int nSlices = sz[dim];
        slices.resize(nSlices);
        region.SetSize(dim, 0);
        
        for (int i = idx[dim]; i < idx[dim] + nSlices; i++) {
            region.SetIndex(dim, i);
            
            FilterType::Pointer filter = FilterType::New();
            filter->SetInput(volume);
            filter->SetExtractionRegion(region);
            filter->SetDirectionCollapseToGuess();
            filter->Update();
            slices[i] = filter->GetOutput();
            slices[i]->DisconnectPipeline();
        }

        return slices;
    }

    RealImage2Vector RealImageTools::computeGaussianSmoothing(RealImage2Vector images, double sigma) {
        RealImage2Vector outputImages;
        typedef itk::RecursiveGaussianImageFilter<RealImage2> FilterType;
        for (int i = 0; i < images.size(); i++) {
            FilterType::Pointer filter = FilterType::New();
            filter->SetInput(images[i]);
            filter->SetSigma(sigma);
            filter->Update();
            outputImages.push_back(filter->GetOutput());
        }
        return outputImages;
    }

    GradientImage2Vector RealImageTools::computeGaussianGradient(RealImage2Vector images, double sigma) {
        GradientImage2Vector outputImages;
        for (int i = 0; i < images.size(); i++) {
            Gradient2FilterType::Pointer filter = Gradient2FilterType::New();
            filter->SetInput(images[i]);
            filter->SetSigma(sigma);
            filter->Update();
            outputImages.push_back(filter->GetOutput());
        }
        return outputImages;
    }
}