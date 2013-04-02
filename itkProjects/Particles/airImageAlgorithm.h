//
//  airImageAlgorithm.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/1/13.
//
//

#ifndef __ParticleGuidedRegistration__airImageAlgorithm__
#define __ParticleGuidedRegistration__airImageAlgorithm__

#include <iostream>
#include <itkImage.h>
#include <itkLabelGeometryImageFilter.h>
#include <itkIsolatedConnectedImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkRegionOfInterestImageFilter.h>
#include <itkMaskImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>


#ifndef NDEBUG
#define __d(msg) std::cout << msg << std::flush
#define __f(msg) std::cout << msg << std::endl
#else
#define __d(msg)
#define __f(msg)
#endif


namespace air {
    template <class T, class L>
    class ImageAlgorithm {
    public:
        typedef typename T::Pointer TPointer;
        typedef typename L::Pointer LPointer;
        typedef typename L::PixelType LPixel;
        typedef typename T::PixelType TPixel;
        typedef typename L::RegionType RegionType;
        typedef typename L::IndexType IndexType;

        typedef itk::ImageRegionIteratorWithIndex<L> LabelIterator;
        typedef itk::ImageRegionIteratorWithIndex<T> ImageIterator;
        typedef itk::LabelGeometryImageFilter<L,T> LabelGeometryFilter;
        typedef itk::RegionOfInterestImageFilter<T,T> ROIFilter;
        typedef itk::MaskImageFilter<T,L,T> MaskFilter;
        typedef itk::BinaryThresholdImageFilter<L,L> BinaryThresholdFilter;

    private:
        typename LabelGeometryFilter::Pointer _labelFilter;

        LPointer ThresholdToBinary(LPointer img) {
            __d("BinaryThreshold() ");
            typename BinaryThresholdFilter::Pointer binThreshFilter = BinaryThresholdFilter::New();
            binThreshFilter->SetInput(img);
            binThreshFilter->SetInsideValue(1);
            binThreshFilter->SetOutsideValue(0);
            binThreshFilter->SetLowerThreshold(1);
            binThreshFilter->SetUpperThreshold(255);
            binThreshFilter->Update();
            __f("done.");
            return binThreshFilter->GetOutput();
        }
        
        TPointer MaskImage(TPointer image, LPointer mask) {
            __d("MaskImage() ");
            typename MaskFilter::Pointer filter = MaskFilter::New();
            filter->SetInput(image);
            filter->SetMaskImage(mask);
            __f("done.");
            return filter->GetOutput();
        }

    public:
        typename LabelGeometryFilter::Pointer GetLabelGeometry(LPointer label, TPointer image) {
            if (label.IsNull()) {
                return typename LabelGeometryFilter::Pointer();
            }
            if (_labelFilter.IsNull() || _labelFilter->GetMTime() < label->GetMTime()) {
                _labelFilter = LabelGeometryFilter::New();
                _labelFilter->SetInput(label);
                _labelFilter->SetIntensityInput(image);
                __d("LabelGeometry() ");
                _labelFilter->Update();
                __f("done.");
            }
            return _labelFilter;
        };

        LPointer ExecuteIsolatedConnectedImageFilter(LPointer labelImg, TPointer srcImg, LPixel fgId, TPixel bgId) {
            typename LabelGeometryFilter::Pointer labelFilter = GetLabelGeometry(labelImg, srcImg);

            RegionType requestedRegion, seed1Region, seed2Region;
            std::vector<LPixel> labels = labelFilter->GetLabels();

            for (int i = 1; i < labels.size(); i++) {
                if (fgId == labels[i]) {
                    seed1Region = labelFilter->GetRegion(fgId);
                } else if (bgId == labels[i]) {
                    seed2Region = labelFilter->GetRegion(bgId);
                }

                // region for whole labels
                if (i == 1) {
                    requestedRegion = labelFilter->GetRegion(labels[i]);
                } else {
                    requestedRegion = UnionRegion(requestedRegion, labelFilter->GetRegion(labels[i]));
                }
            }
            srcImg->SetRequestedRegion(requestedRegion);
            std::cout << "RequestedRegion: " << requestedRegion << std::endl;
            std::cout << "Seed1 Region: " << seed1Region << std::endl;
            std::cout << "Seed2 Region: " << seed2Region << std::endl;

//
//            __d("ROI Filter running");
//            typename ROIFilter::Pointer roiFilter = ROIFilter::New();
//            roiFilter->SetInput(srcImg);
//            roiFilter->SetRegionOfInterest(requestedRegion);
//            __f(" done.");

//            TPointer roiImage = roiFilter->GetOutput();
//            roiImage->Print(std::cout);

            TPointer roiImage = srcImg;
//            roiImage = MaskImage(srcImg, labelImg);

            typedef itk::IsolatedConnectedImageFilter<T,L> FilterType;
            typename FilterType::Pointer filter = FilterType::New();
            filter->SetInput(roiImage);
            filter->SetReplaceValue(fgId);


            // add seed1
            LabelIterator seed1Iter(labelImg, seed1Region);
            for (seed1Iter.GoToBegin(); !seed1Iter.IsAtEnd(); ++seed1Iter) {
                if (seed1Iter.Get() == fgId) {
                    filter->AddSeed1(seed1Iter.GetIndex());
                }
            }

            // add seed2
            LabelIterator seed2Iter(labelImg, seed2Region);
            for (seed2Iter.GoToBegin(); !seed2Iter.IsAtEnd(); ++seed2Iter) {
                if (seed2Iter.Get() == bgId) {
                    filter->AddSeed2(seed2Iter.GetIndex());
                }
            }

            filter->FindUpperThresholdOff();
            try {
                __d("IsoRG running " << std::flush);
                filter->Update();
                __f(" done.");
            } catch (itk::ExceptionObject& ex) {
                ex.Print(std::cout);
            }

            if (filter->GetThresholdingFailed()) {
                std::cout << "Filter Execution Failed" << std::endl;
                return LPointer();
            }

            LPointer labelOut = filter->GetOutput();
            labelOut->DisconnectPipeline();

            // roll back original image
            srcImg->SetRequestedRegion(srcImg->GetBufferedRegion());

            return labelOut;
        }

        RegionType UnionRegion(RegionType a, RegionType b) {
            IndexType lowerA = a.GetIndex();
            IndexType upperA = a.GetUpperIndex();
            IndexType lowerB = b.GetIndex();
            IndexType upperB = b.GetUpperIndex();
            IndexType lowerOut, upperOut;
            for (int i = 0; i < IndexType::GetIndexDimension(); i++) {
                lowerOut[i] = std::min(lowerA[i], lowerB[i]);
                upperOut[i] = std::max(upperA[i], upperB[i]);
            }
            RegionType region;
            region.SetIndex(lowerOut);
            region.SetUpperIndex(upperOut);
            return region;
        }

        
    };
}

#endif /* defined(__ParticleGuidedRegistration__airImageAlgorithm__) */
