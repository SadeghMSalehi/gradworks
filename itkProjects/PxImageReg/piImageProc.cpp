//
//  piImageProc.cpp
//  PxImageReg
//
//  Created by Joohwi Lee on 11/5/13.
//
//

#include "piImageProc.h"
#include "piImageIO.h"

#include <itkSignedDanielssonDistanceMapImageFilter.h>
#include <itkExtractImageFilter.h>
#include <itkVectorMagnitudeImageFilter.h>

using namespace std;

namespace pi {

    typedef itk::SignedDanielssonDistanceMapImageFilter<LabelImage, RealImage> SignedDistanceMapFilterType;

    class OffsetToVector {
    public:
        bool operator!=(const OffsetToVector& o) {
            return false;
        }
        bool operator==(const OffsetToVector& o) {
            return true;
        }
        inline VectorType operator()(SignedDistanceMapFilterType::VectorImageType::PixelType o) {
            VectorType v;
            fordim (k) {
                v[k] = o[k];
            }
            return v;
        }
    };

    VectorImage::Pointer ComputeDistanceMap(LabelImage::Pointer img, std::string magnitudeFile) {
        cout << "Computing distance map ..." << flush;
        LabelImage::Pointer binaryMap = img;

        // construct signed distance filter
        SignedDistanceMapFilterType::Pointer distmapFilter = SignedDistanceMapFilterType::New();
        distmapFilter->SetInput(binaryMap);
        distmapFilter->InsideIsPositiveOff();
        distmapFilter->UseImageSpacingOn();
        distmapFilter->Update();
        SignedDistanceMapFilterType::OutputImagePointer distmap = distmapFilter->GetDistanceMap();

        if (magnitudeFile != "") {
            ImageIO<SignedDistanceMapFilterType::OutputImageType> io;
            io.WriteImage(magnitudeFile, distmap);
        }


        /** Pointer Type for the vector distance image */
        OffsetImage::Pointer distanceOffsetImage = distmapFilter->GetVectorDistanceMap();
        ImageIO<VectorImage> io;
        VectorImage::Pointer distanceVectorImage = io.NewImageS<LabelImage>(binaryMap);

        const int nPixels = img->GetPixelContainer()->Size();

        OffsetType* offsetPointer = distanceOffsetImage->GetBufferPointer();
        VectorType* vectorPointer = distanceVectorImage->GetBufferPointer();

        for (int i = 0; i < nPixels; i++) {
            fordim (k) {
                vectorPointer[i][k] = offsetPointer[i][k];
            }
        }
        return distanceVectorImage;
    }



    GradientImage::Pointer ComputeGaussianGradient(LabelImage::Pointer img, double sigma) {
        typedef itk::CastImageFilter<LabelImage ,RealImage > CastToRealFilterType;

        CastToRealFilterType::Pointer toReal = CastToRealFilterType::New();
        toReal->SetInput(img);
        toReal->Update();

        GaussianGradientFilterType::Pointer gradientFilter = GaussianGradientFilterType::New();
        gradientFilter->SetInput(toReal->GetOutput());
        if (sigma > -1) {
            gradientFilter->SetSigma(sigma);
        }
        gradientFilter->Update();
        return gradientFilter->GetOutput();
    }


    GradientImage::Pointer ComputeGaussianGradient(RealImage::Pointer img, double sigma) {
        GaussianGradientFilterType::Pointer gradientFilter = GaussianGradientFilterType::New();
        gradientFilter->SetInput(img);
        if (sigma > -1) {
            gradientFilter->SetSigma(sigma);
        }
        gradientFilter->Update();
        return gradientFilter->GetOutput();
    }


    RealImage::Pointer ComputeGradientMagnitude(GradientImage::Pointer image) {
        typedef itk::VectorMagnitudeImageFilter<GradientImage, RealImage> MagnitudeFilterType;
        MagnitudeFilterType::Pointer magFilter = MagnitudeFilterType::New();
        magFilter->SetInput(image);
        magFilter->Update();
        return magFilter->GetOutput();
    }

    LabelImage3::Pointer CreateImage3(LabelImage::Pointer refImage, int m) {
        ImageIO<LabelImage3> io;
        LabelImage::SizeType sz = refImage->GetBufferedRegion().GetSize();
        LabelImage3::Pointer tracker = io.NewImageT(sz[0], sz[1], m);
        LabelImage3::SpacingType spacing;
        LabelImage3::PointType origin;
        LabelImage3::DirectionType direction;
        direction.Fill(0);
        fordim (k) {
            spacing[k] = refImage->GetSpacing()[k];
            origin[k] = refImage->GetOrigin()[k];
            fordim (l) {
                direction[k][l] = refImage->GetDirection()[k][l];
            }
        }
        spacing[2] = spacing[0];
        origin[2] = 0;
        direction[2][2] = 1;

        tracker->FillBuffer(0);
        tracker->SetSpacing(spacing);
        tracker->SetOrigin(origin);
        tracker->SetDirection(direction);

        return tracker;
    }
    
    RealImage2Vector SliceVolume(RealImage3::Pointer volume, int dim) {
        RealImage2Vector slices;
        RealImage3::RegionType region = volume->GetBufferedRegion();
        RealImage3::SizeType sz = region.GetSize();
        const RealImage3::IndexType idx = region.GetIndex();
        
        const int nSlices = sz[dim];
        slices.resize(nSlices);
        region.SetSize(dim, 0);
        
        typedef itk::ExtractImageFilter<RealImage3, RealImage2> FilterType;
        
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
}