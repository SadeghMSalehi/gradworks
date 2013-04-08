//
//  piContourSystem.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/6/13.
//
//

#include <itkResampleImageFilter.h>

#include "piImageIO.h"
#include "piContourSystem.h"
#include "piEntropyComputer.h"

namespace pi {
    void ContourSystem::AllocateAttributeBuffer() {
        int nElems = (2 * _width + 1) * (2 * _height + 1);

        _attrs = new EntropyComputer<float>(2, _particles.size(), nElems);
        _attrs->dataIter.FirstData();
    }


    void ContourSystem::ExtractAttributes() {
        _attrs->dataIter.FirstData();
        for (int i = 0; i < _sliceImages.size(); i++) {
            for (int j = 0; j < _particles.size(); j++) {
                float* attrs = _attrs->dataIter.Sample();
                ExtractParticleAttributes(_sliceImages[i], _particles[j], attrs, _width, _height);
                _attrs->dataIter.NextSample();
            }
            _attrs->dataIter.NextData();
        }
    }

    // Extract neighborhood intensity values among particles
    void ContourSystem::ExtractParticleAttributes(RealSlice::Pointer sliceImg,
                                                  Particle& par, float* attrOut, int nWidth, int nHeight) {
        RealSlice::PixelContainer::Pointer importer = RealSlice::PixelContainer::New();
        importer->SetContainerManageMemory(false);

        // set poitner for i-th sample
        importer->SetImportPointer(attrOut, (2*nWidth+1) * (2*nHeight+1));

        RealSlice::RegionType region;
        region.SetSize(0, nWidth);
        region.SetSize(1, nHeight);

        RealSlice::SpacingType spacing = sliceImg->GetSpacing();
        RealSlice::Pointer attrImage = RealSlice::New();
        attrImage->SetPixelContainer(importer);
        attrImage->SetSpacing(spacing);

        attrImage->SetRequestedRegion(region);
        attrImage->SetBufferedRegion(region);
        attrImage->SetLargestPossibleRegion(region);

        typedef itk::AffineTransform<double,2> TransformType;
        TransformType::Pointer transform = TransformType::New();
        TransformType::OutputVectorType offset;
        offset[0] = (par.x[0] - nWidth) * spacing[0];
        offset[1] = (par.x[1] - nHeight) * spacing[1];

        transform->SetIdentity();
        transform->Translate(offset);

        typedef itk::ResampleImageFilter<RealSlice, RealSlice> Resampler;
        Resampler::Pointer resampler = Resampler::New();
        resampler->SetInput(sliceImg);
        resampler->SetTransform(transform.GetPointer());
        resampler->SetReferenceImage(attrImage);
        resampler->UseReferenceImageOn();
        resampler->GraftOutput(attrImage);
        resampler->Update();
        RealSlice::Pointer output = resampler->GetOutput();
    }

    void ContourSystem::Track(int nSlices) {
        for (int i = 1; i <= nSlices && _sliceImages.size(); i++) {
            _attrs->dataIter.At(_sliceIdx + i, 0);
            for (int j = 0; j < _particles.size(); j++) {
                float* attrs = _attrs->dataIter.Sample();
                ExtractParticleAttributes(_sliceImages[_sliceIdx + i], _particles[j], attrs, _width, _height);
                _attrs->dataIter.NextSample();
            }
            _attrs->MoveToCenter();
            _attrs->ComputeCovariance(1);
            _attrs->ComputeGradient();
        }
    }
}