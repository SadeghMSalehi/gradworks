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

#include "itkMyScaleVersor3DTransformOptimizer.h"
#include "itkScaleVersor3DTransform.h"
#include "itkMyRegistration.h"
#include "QThread"

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

    ScaleRegistration::Pointer _registrationAlgorithm;

    int CurrentSliceIndex;

    itkMyCore() {
        _minSliceIndex = 0;
        _maxSliceIndex = 0;
        CurrentSliceIndex = 0;
    }

    int GetMinSliceIndex() {
        return _minSliceIndex;
    }

    int GetMaxSliceIndex() {
        return _maxSliceIndex;
    }

    void LoadImage(const char* file) {
        try {
            cout << "LoadImage() >> Current Thread Id: " << QThread::currentThreadId() << endl;

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
            InverseLabelSlice = LabelSliceType::New();
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

    }

    void ApplyTransform(int historyId) {
        if (_registrationAlgorithm.IsNotNull()) {
            InverseLabelSlice->SetLabel(_registrationAlgorithm->TransformFixedLabel(historyId));
        }
    }

    void WriteLastTransform(const char* fileName) {
        _registrationAlgorithm->WriteTransform(fileName, -1);
    }

    void ApplyLastTransform() {
        int historyId = _registrationAlgorithm->GetTransformHistory().size() - 1;
        if (historyId < 0) {
            return;
        }
        ApplyTransform(historyId);
    }

    void PrepareRegistration() {
        _registrationAlgorithm = ScaleRegistration::New();
    }
    
    ScaleRegistration::TransformHistoryType RunRegistration() {
        if (_registrationAlgorithm.IsNull()) {
            return ScaleRegistration::TransformHistoryType();
        }
        _registrationAlgorithm->SetImages(SourceSlice->GetSourceImage(), LabelSlice->GetViewImage(), TargetSlice->GetSourceImage());
        try {
            _registrationAlgorithm->RunRegistration();
        } catch (...) {
            //cout << "I ate the registratione exception! [" << QThread::currentThreadId() << "]; RegistrationAlgorithm " << _registrationAlgorithm.GetPointer() <<  endl;
        }
        return _registrationAlgorithm->GetTransformHistory();
    }
};

#endif /* defined(__itktools__itkMyCore__) */