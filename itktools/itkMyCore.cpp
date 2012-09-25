//
//  itkMyCore.cpp
//  itktools
//
//  Created by Joohwi Lee on 9/20/12.
//
//

#include "itkMyCore.h"
#include "QVector"

int itkMyCore::GetMinSliceIndex() {
	return _minSliceIndex;
}

int itkMyCore::GetMaxSliceIndex() {
	return _maxSliceIndex;
}

void itkMyCore::LoadImage(const char* file) {
	try {
		cout << "LoadImage() >> Current Thread Id: "
				<< QThread::currentThreadId() << endl;

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

void itkMyCore::LoadTarget(const char* file) {
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

void itkMyCore::LoadLabelIfGrayImageLoaded(const char* file) {
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

void itkMyCore::SetCurrentSlice(int slice) {
	if (slice < _minSliceIndex || slice >= _maxSliceIndex) {
		CurrentSliceIndex = SourceSlice->ComputeSliceAtCenter();
	} else {
		CurrentSliceIndex = slice;
	}
}

void itkMyCore::LoadTransform(const char* fileName) {

}

void itkMyCore::ApplyTransform(int historyId) {
	if (_registrationAlgorithm.IsNotNull()) {
		InverseLabelSlice->SetLabel(
				_registrationAlgorithm->TransformFixedLabel(historyId));
	}
}

void itkMyCore::WriteLastTransform(const char* fileName) {
	_registrationAlgorithm->WriteTransform(fileName, -1);
}

void itkMyCore::ApplyLastTransform() {
	int historyId = _registrationAlgorithm->GetTransformHistory().size() - 1;
	if (historyId < 0) {
		return;
	}
	ApplyTransform(historyId);
}

void itkMyCore::PrepareRegistration() {
	_registrationAlgorithm = RegistrationMethod::New();
}

RegistrationMethod::TransformHistoryType itkMyCore::RunRegistration() {
	if (_registrationAlgorithm.IsNull()) {
		return RegistrationMethod::TransformHistoryType();
	}
	_registrationAlgorithm->SetImages(SourceSlice->GetSourceImage(),
			LabelSlice->GetViewImage(), TargetSlice->GetSourceImage());
	try {
		_registrationAlgorithm->RunRegistration();
	} catch (...) {
		//cout << "I ate the registratione exception! [" << QThread::currentThreadId() << "]; RegistrationAlgorithm " << _registrationAlgorithm.GetPointer() <<  endl;
	}
	return _registrationAlgorithm->GetTransformHistory();
}
