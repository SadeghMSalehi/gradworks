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
//#include "itkMySlicer.h"
#include "itkMySlice.h"
#include "itkImageIO.h"
#include "itkExceptionObject.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkTransformFileReader.h"
#include "itkResampleImageFilter.h"

#include "itkMyScaleVersor3DTransformOptimizer.h"
#include "itkScaleVersor3DTransform.h"
#include "itkMultiRegMethod.h"
//#include "itkMyRegistration.h"
#include "QThread"
#include "QImage"
//

class itkMyCore {
private:
    BitmapType::Pointer _sourceBitmap[ImageType::ImageDimension];
    BitmapType::Pointer _targetBitmap[ImageType::ImageDimension];
    BitmapType::Pointer _labelBitmap[ImageType::ImageDimension];
    BitmapType::Pointer _transformedLabelBitmap[ImageType::ImageDimension];

    ImageSlicer::Pointer _sourceSlicer;
    ImageSlicer::Pointer _targetSlicer;
    LabelSlicer::Pointer _labelSlicer;
    LabelSlicer::Pointer _transformedLabelSlicer;

    int _labelOpacity;
    void UpdateBitmap(ImageKindType kind, int dir, int opacity);
    
public:
	itkcmds::itkImageIO<ImageType> imageIO;
	itkcmds::itkImageIO<LabelType> labelIO;

	RegistrationMethod::Pointer _registrationAlgorithm;

	itkMyCore() {
        _labelOpacity = 120;
	}

//	int GetMinSliceIndex();
//	int GetMaxSliceIndex(int xyz = 2);
//    int GetCurrentSliceIndex(int xyz);

	void LoadImage(const char* file);
	void LoadTarget(const char* file);
	void LoadLabel(const char* file);

    void SetCurrentSlice(int dir, int sliceIdx);

	ImageType::IndexType GetCurrentSliceIndex() {
		return _sourceSlicer->GetCurrentSliceIndex();
	}

    ImageType::SizeType GetMaxSliceIndex() {
        return _sourceSlicer->GetMaxSliceIndex();
    }

	void PrepareRegistration();
	RegistrationMethod::TransformHistoryType RunRegistration();
	void ApplyLastTransform();
	void ApplyTransform(int historyId);
	void WriteLastTransform(const char* fileName);

    BitmapType::Pointer GetSourceBitmap(int dir);
    BitmapType::Pointer GetTargetBitmap(int dir);
    BitmapType::Pointer GetLabelBitmap(int dir);
    BitmapType::Pointer GetTransformedBitmap(int dir);

    QImage ConvertToQImage(ImageKindType kind, int dir);
    
    void SetLabelOpacity(int opacity) {
        _labelOpacity = opacity;
    }

	void LoadTransform(const char* fileName) {
        this->PrepareRegistration();
        _registrationAlgorithm->LoadParameterHistory(fileName);
    }
    
    void SaveTransform(const char* filename) {
        if (_registrationAlgorithm.IsNotNull()) {
            _registrationAlgorithm->SaveParameterHistory(filename);
        }
    }


    void ExecuteCommandLine(int argc, char* argv[]);
};

#endif /* defined(__itktools__itkMyCore__) */
