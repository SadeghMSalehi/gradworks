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

    VectorImage::Pointer ComputeDistanceMap(LabelImage::Pointer img) {
        cout << "Computing distance map ..." << flush;
        LabelImage::Pointer binaryMap = img;

        // construct signed distance filter
        SignedDistanceMapFilterType::Pointer distmapFilter = SignedDistanceMapFilterType::New();
        distmapFilter->SetInput(binaryMap);
        distmapFilter->InsideIsPositiveOff();
        distmapFilter->UseImageSpacingOn();
        distmapFilter->Update();
        SignedDistanceMapFilterType::OutputImagePointer distmap = distmapFilter->GetDistanceMap();

        typedef
        itk::UnaryFunctorImageFilter<SignedDistanceMapFilterType::VectorImageType, VectorImage, OffsetToVector> OffsetToVectorCastFilterType;

        OffsetToVectorCastFilterType::Pointer caster = OffsetToVectorCastFilterType::New();
        caster->SetInput(distmapFilter->GetVectorDistanceMap());
        caster->Update();
        cout << "done" << endl;
        return caster->GetOutput();
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
}