//
//  itkMyCore.h
//  itktools
//
//  Created by Joohwi Lee on 9/20/12.
//
//

#ifndef __itktools__itkMyCore__
#define __itktools__itkMyCore__

#include <iostream>
#include "itkImageIO.h"
#include "itkExceptionObject.h"
#include "itkExtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

typedef itk::Image<unsigned short, 3> ImageType;
typedef itk::Image<unsigned short, 2> SliceType;
typedef itk::Image<int, 2> BitmapType;

class itkMyCore {
public:
    const char* ImageName;

    itkcmds::itkImageIO<ImageType> imageIO;
    ImageType::Pointer GrayImage;
    ImageType::SizeType GrayImageSize;
    BitmapType::Pointer CurrentSlice;

    void LoadImage(const char* file) {
        try {
            ImageName = file;
            GrayImage = imageIO.ReadImageT(ImageName);
            GrayImageSize = GrayImage->GetBufferedRegion().GetSize();
            ExtractSlice();
        } catch (itk::ExceptionObject& ex) {
            cout << ex << endl;
        }
    }

    void ExtractSlice(int slice = -1) {
        typedef itk::ExtractImageFilter<ImageType, SliceType> SliceExtractor;
        typedef itk::RescaleIntensityImageFilter<SliceType, BitmapType> IntensityFilter;
        SliceExtractor::Pointer slicer = SliceExtractor::New();
        ImageType::RegionType sliceRegion;
        ImageType::SizeType sliceSize;
        ImageType::IndexType sliceIndex;
        sliceSize[0] = GrayImageSize[0];
        sliceSize[1] = GrayImageSize[1];
        sliceSize[2] = 0;
        sliceIndex[0] = sliceIndex[1] = 0;
        sliceIndex[2] = (slice < 0 || slice >= GrayImageSize[2]) ? (GrayImageSize[2] / 2) : slice;
        sliceRegion.SetSize(sliceSize);
        sliceRegion.SetIndex(sliceIndex);
        slicer->SetInput(GrayImage);
        slicer->SetExtractionRegion(sliceRegion);
        slicer->SetDirectionCollapseToSubmatrix();
        slicer->Update();
        SliceType::Pointer tmp = slicer->GetOutput();
        IntensityFilter::Pointer filter = IntensityFilter::New();
        filter->SetInput(tmp);
        filter->SetOutputMinimum(0);
        filter->SetOutputMaximum(255);
        filter->Update();
        CurrentSlice = filter->GetOutput();
    }
};

#endif /* defined(__itktools__itkMyCore__) */
