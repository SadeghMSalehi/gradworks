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
#include "itkMyTypes.h"
#include "itkMySlicer.h"
#include "itkImageIO.h"
#include "itkExceptionObject.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkTransformFileReader.h"
#include "itkResampleImageFilter.h"

class itkMyCore {
private:
    int _maxSliceIndex;
    int _minSliceIndex;

public:

    itkcmds::itkImageIO<ImageType> imageIO;
    itkcmds::itkImageIO<LabelType> labelIO;

    GraySliceType::Pointer SourceSlice;
    GraySliceType::Pointer TargetSlice;
    LabelSliceType::Pointer LabelSlice;
    LabelSliceType::Pointer InverseLabelSlice;

    TransformReaderType::TransformListType* TransformList;

    int CurrentSliceIndex;

    itkMyCore() {
        _minSliceIndex = 0;
        _maxSliceIndex = 0;
        CurrentSliceIndex = 0;
        TransformList = NULL;
    }

    int GetMinSliceIndex() {
        return _minSliceIndex;
    }

    int GetMaxSliceIndex() {
        return _maxSliceIndex;
    }

    void LoadImage(const char* file) {
        try {
            ImageType::Pointer image = imageIO.ReadImageT(file);
            SourceSlice = GraySliceType::New();
            SourceSlice->SetName(std::string(file));
            SourceSlice->SetInput(image);
            CurrentSliceIndex = SourceSlice->ComputeSliceAtCenter();

            _minSliceIndex = 0;
            _maxSliceIndex = SourceSlice->GetSize()[2];
        } catch (itk::ExceptionObject& ex) {
            cout << ex << endl;
        }
    }

    void LoadTarget(const char* file) {
        try {
            ImageType::Pointer image = imageIO.ReadImageT(file);
            TargetSlice = GraySliceType::New();
            TargetSlice->SetName(std::string(file));
            TargetSlice->SetInput(image);
            TargetSlice->UpdateSlice(CurrentSliceIndex, 0xff);
        } catch (itk::ExceptionObject& ex) {
            cout << ex << endl;
        }
    }

    void LoadLabelIfGrayImageLoaded(const char* file) {
        try {
            LabelType::Pointer image = labelIO.ReadImageT(file);
            LabelSlice = LabelSliceType::New();
            LabelSlice->SetName(std::string(file));
            LabelSlice->SetLabel(image);
            LabelSlice->UpdateSlice(CurrentSliceIndex, 100);
        } catch (itk::ExceptionObject& ex) {
            cout << ex << endl;
        }
    }

    void SetCurrentSlice(int slice) {
        if (slice < _minSliceIndex || slice >= _maxSliceIndex) {
            CurrentSliceIndex = SourceSlice->ComputeSliceAtCenter();
        } else {
            CurrentSliceIndex = slice;
        }
    }

    void LoadTransform(const char* fileName) {
        itk::TransformFileReader::Pointer reader = itk::TransformFileReader::New();
        reader->SetFileName(fileName);
        reader->Update();
        TransformList = reader->GetTransformList();
    }

    void ApplyTransform() {
        if (TransformList == NULL || LabelSlice->GetViewImage().IsNull()) {
            return;
        }
        TransformBase::Pointer transform = TransformList->front();
        MatrixTransformType* matrixTransform = dynamic_cast<MatrixTransformType*>(transform.GetPointer());
        if (matrixTransform == NULL) {
            return;
        }
        MatrixTransformType::InverseTransformBasePointer inverseTransform = matrixTransform->GetInverseTransform();

        typedef itk::ResampleImageFilter<LabelType, LabelType> ResampleFilter;
        ResampleFilter::Pointer resampler = ResampleFilter::New();
        resampler->SetInput(LabelSlice->GetViewImage());
        resampler->SetTransform(inverseTransform);
        resampler->SetInterpolator(InterpolatorNN::New());
        resampler->SetReferenceImage(TargetSlice->GetViewImage());
        resampler->SetUseReferenceImage(true);
        resampler->Update();
        InverseLabelSlice = LabelSliceType::New();
        InverseLabelSlice->SetLabel(resampler->GetOutput());
    }
};

#endif /* defined(__itktools__itkMyCore__) */