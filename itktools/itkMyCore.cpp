//
//  itkMythis->cpp
//  itktools
//
//  Created by Joohwi Lee on 9/20/12.
//
//

#include "itkMyCore.h"
#include "itkMySlice.h"
#include "QVector"

QImage itkMyCore::ConvertToQImage(ImageKindType kind, int dir) {
    unsigned char* buffer = NULL;
    if (kind == Source) {
        buffer = (unsigned char*) GetSourceBitmap(dir)->GetBufferPointer();
    } else if (kind == Target) {
        buffer = (unsigned char*) GetTargetBitmap(dir)->GetBufferPointer();
    } else if (kind == Label) {
        buffer = (unsigned char*) GetLabelBitmap(dir)->GetBufferPointer();
    } else if (kind == TransformedLabel) {
        buffer = (unsigned char*) GetTransformedBitmap(dir)->GetBufferPointer();
    }
    return QImage((unsigned char*) buffer, _sourceSlicer->GetSize(dir)[0], _sourceSlicer->GetSize(dir)[1], QImage::Format_ARGB32);
}

void itkMyCore::UpdateBitmap(ImageKindType kind, int dir, int opacity) {
    if (kind == 0) {
        GrayImageFilter::Pointer filter = GrayImageFilter::New();
        filter->SetInput(_sourceSlicer->GetSlicePointer(dir));
        filter->GetFunctor().SetImageStatistics(_sourceSlicer->GetMinIntensity(), _sourceSlicer->GetMaxIntensity(), _sourceSlicer->GetMeanIntensity(), _sourceSlicer->GetStdIntensity());
        filter->GetFunctor().SetInsideAlpha(opacity);
        filter->Update();
        _sourceBitmap[dir] = filter->GetOutput();
    } else if (kind == 1) {
        GrayImageFilter::Pointer filter = GrayImageFilter::New();
        filter->SetInput(_targetSlicer->GetSlicePointer(dir));
        filter->GetFunctor().SetImageStatistics(_targetSlicer->GetMinIntensity(), _targetSlicer->GetMaxIntensity(), _targetSlicer->GetMeanIntensity(), _targetSlicer->GetStdIntensity());
        filter->GetFunctor().SetInsideAlpha(opacity);        
        filter->Update();
        _targetBitmap[dir] = filter->GetOutput();
    } else if (kind == Label){
        LabelMapFilter::Pointer filter = LabelMapFilter::New();
        filter->SetInput(_labelSlicer->GetSlicePointer(dir));
        filter->GetFunctor().SetInsideAlpha(opacity);
        filter->Update();
        _labelBitmap[dir] = filter->GetOutput();
    } else if (kind == TransformedLabel) {
        if (_transformedLabelSlicer.IsNull()) {
            return;
        }
        LabelMapFilter::Pointer filter = LabelMapFilter::New();
        filter->SetInput(_transformedLabelSlicer->GetSlicePointer(dir));
        filter->GetFunctor().SetInsideAlpha(opacity);
        filter->Update();
        _transformedLabelBitmap[dir] = filter->GetOutput();
    }
}

BitmapType::Pointer itkMyCore::GetSourceBitmap(int dir) {
    if (_sourceBitmap[dir].IsNull() || _sourceSlicer->GetMTime() > _sourceBitmap[dir]->GetMTime()) {
        UpdateBitmap(Source, dir, 255);
    }
    return _sourceBitmap[dir];
}

BitmapType::Pointer itkMyCore::GetTargetBitmap(int dir) {
    if (_targetBitmap[dir].IsNull() || _targetSlicer->GetMTime() > _targetBitmap[dir]->GetMTime()) {
        UpdateBitmap(Target, dir, 255);
    }
    return _targetBitmap[dir];
}

BitmapType::Pointer itkMyCore::GetLabelBitmap(int dir) {
    UpdateBitmap(Label, dir, _labelOpacity);
    return _labelBitmap[dir];
}

BitmapType::Pointer itkMyCore::GetTransformedBitmap(int dir) {
    UpdateBitmap(TransformedLabel, dir, _labelOpacity);
    return _transformedLabelBitmap[dir];
}

void itkMyCore::LoadImage(const char* file) {
	try {
		ImageType::Pointer image = imageIO.ReadImageT(file);
        _sourceSlicer = ImageSlicer::New();
        _sourceSlicer->SetVolume(image);
	} catch (itk::ExceptionObject& ex) {
		cout << ex << endl;
	}
}

void itkMyCore::LoadTarget(const char* file) {
	try {
		ImageType::Pointer image = imageIO.ReadImageT(file);
        if (image->GetBufferedRegion().GetSize() != _sourceSlicer->GetMaxSliceIndex()) {
            image = imageIO.ResampleImageAs(image, _sourceSlicer->GetVolume());
        }
        _targetSlicer = ImageSlicer::New();
        _targetSlicer->SetVolume(image);
	} catch (itk::ExceptionObject& ex) {
		cout << ex << endl;
	}
}

void itkMyCore::LoadLabel(const char* file) {
	try {
		LabelType::Pointer image = labelIO.ReadImageT(file);
        _labelSlicer = LabelSlicer::New();
        _labelSlicer->SetVolume(image);
	} catch (itk::ExceptionObject& ex) {
		cout << ex << endl;
	}
}

void itkMyCore::SetCurrentSlice(int dir, int sliceIdx) {
    if (_sourceSlicer.IsNotNull()) {
        _sourceSlicer->Update(dir, sliceIdx);
    }
    if (_targetSlicer.IsNotNull()) {
        _targetSlicer->Update(dir, sliceIdx);
    }
    if (_labelSlicer.IsNotNull()) {
        _labelSlicer->Update(dir, sliceIdx);
    }
}


void itkMyCore::ApplyTransform(int historyId) {
    int historyCount = _registrationAlgorithm->GetTransformHistory().size();
	if (_registrationAlgorithm.IsNotNull()) {
        if (historyId < 0) {
            historyId += _registrationAlgorithm->GetTransformHistory().size();
        } else if (historyId >= historyCount) {
            historyId = historyCount - 1;
        }
        _transformedLabelSlicer = LabelSlicer::New();
        _transformedLabelSlicer->SetVolume(_registrationAlgorithm->TransformFixedLabel(historyId));
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
	_registrationAlgorithm->SetImages(_sourceSlicer->GetVolume(),
                                      _labelSlicer->GetVolume(), _targetSlicer->GetVolume());
}

RegistrationMethod::TransformHistoryType itkMyCore::RunRegistration() {
	if (_registrationAlgorithm.IsNull()) {
		return RegistrationMethod::TransformHistoryType();
	}
	try {
		_registrationAlgorithm->RunRegistration();
	} catch (...) {
		cout << "I ate a registratione exception! [" << QThread::currentThreadId() << "]; RegistrationAlgorithm " << _registrationAlgorithm.GetPointer() <<  endl;
	}
	return _registrationAlgorithm->GetTransformHistory();
}

void itkMyCore::ExecuteCommandLine(int argc, char **argv) {
	string cmd(argv[1]);
	if (cmd == "registration") {
		this->LoadImage(argv[2]);
		this->LoadTarget(argv[3]);
		this->LoadLabel(argv[4]);
		this->PrepareRegistration();
		this->RunRegistration();
		this->ApplyLastTransform();
//		if (this->InverseLabelSlice.IsNotNull()) {
//			LabelType::Pointer label = this->InverseLabelSlice->GetViewImage();
//			if (label.IsNotNull()) {
//				itkcmds::itkImageIO<LabelType> io;
//				io.WriteImageT(argv[5], label);
//			}
//		}
		if (argc > 6) {
			this->WriteLastTransform(argv[6]);
		}
	} else if (cmd == "slice") {
        this->LoadImage(argv[2]);

    } else if (cmd == "test") {
    }

}

