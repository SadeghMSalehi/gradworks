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
#include "itkMultiRegMethod.h"
//#include "itkMyRegistration.h"
#include "QThread"

//

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

	RegistrationMethod::Pointer _registrationAlgorithm;

	int CurrentSliceIndex;

	itkMyCore() {
		_minSliceIndex = 0;
		_maxSliceIndex = 0;
		CurrentSliceIndex = 0;
	}

	int GetMinSliceIndex();
	int GetMaxSliceIndex();

	void LoadImage(const char* file);
	void LoadTarget(const char* file);
	void LoadLabelIfGrayImageLoaded(const char* file);

	void SetCurrentSlice(int slice);
	void LoadTransform(const char* fileName);
	void ApplyTransform(int historyId);
	void WriteLastTransform(const char* fileName);

	void ApplyLastTransform();
	void PrepareRegistration();

	RegistrationMethod::TransformHistoryType RunRegistration();
};

#endif /* defined(__itktools__itkMyCore__) */
