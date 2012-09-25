/*
 * itkMultiRegMethod.cpp
 *
 *  Created on: Sep 24, 2012
 *      Author: joohwi
 */

#include "itkMultiRegMethod.h"
#include "itkContinuousIndex.h"

void itkMultiRegMethod::Execute(itk::Object* caller,
		const itk::EventObject& event) {
	Execute((const Object*) caller, event);
}

void itkMultiRegMethod::OnIterate(const OptimizerType* optimizer) {
	if (++_numOfIterations % _updateInterval == 0) {
		itk::RealTimeClock::TimeStampType t = _clock->GetTimeInSeconds();
		std::cout << _numOfIterations << "\t" << optimizer->GetValue() << "\t"
				<< optimizer->GetCurrentPosition() << "\t" << (t - _lastTime)
				<< " secs" << std::endl;
		_transformHistory.push_back(optimizer->GetCurrentPosition());

		_lastTime = t;
	}
}

void itkMultiRegMethod::Execute(const itk::Object* object,
		const itk::EventObject& event) {
	if (object == NULL) {
		std::cout << "Null sender is not processed..." << std::endl;
		return;
	}
	if (typeid(event) == typeid(itk::IterationEvent)) {
		const OptimizerType* optimizer =
				dynamic_cast<const OptimizerType*>(object);
		if (optimizer == NULL) {
			cout << "Wrong optimizer type" << endl;
			return;
		}
		OnIterate(optimizer);
	} else if (typeid(event) == typeid(itk::StartEvent)) {
		std::cout << "Optimization has started ..." << std::endl;
		_lastTime = _clock->GetTimeInSeconds();
	} else if (typeid(event) == typeid(itk::EndEvent)) {
		std::cout << "Optimization has finisehd ..." << std::endl;
	}
}

void itkMultiRegMethod::Update() {
	this->Execute((const Object*) NULL, itk::IterationEvent());
}

void itkMultiRegMethod::SetImages(ImageType::Pointer source,
		LabelType::Pointer sourceLabel, ImageType::Pointer target) {
	_movingImage = target;
	_fixedImage = source;
	_fixedImageLabel = sourceLabel;
}

void itkMultiRegMethod::RunRegistration() {
	cout << "Generating label indexes ..." << endl;

	_centerOfRotation.resize(30);
	_labelIndexes.resize(30);

	int nPixels[30];
	for (int i = 0; i < 30; i++) {
		nPixels[i] = 0;
		_centerOfRotation[i].SetSize(3);
	}

	itk::ImageRegionConstIteratorWithIndex<LabelType> labelIter(
			_fixedImageLabel, _fixedImageLabel->GetBufferedRegion());

	int maxLabel = 0;
	for (labelIter.GoToBegin(); !labelIter.IsAtEnd(); ++labelIter) {
		maxLabel ++;
		LabelType::PixelType label = labelIter.Get();
		if (label > 0) {
			LabelType::IndexType idx = labelIter.GetIndex();
			_labelIndexes[label].push_back(idx);
			for (int i = 0; i < ImageType::ImageDimension; i++) {
				_centerOfRotation[label][i] += idx[i];
			}
			nPixels[label]++;
		}
	}
	for (int label = 0; label < maxLabel; label++) {
		for (int i = 0; i < ImageType::ImageDimension; i++) {
			_centerOfRotation[label][i] /= nPixels[label];
		}
	}

	cout << "Registration starting ..." << endl;

	typedef itk::MyMetric<ImageType,ImageType> SubMetricType;
	typedef MetaMetrics<SubMetricType> MetaMetricType;

	MetaMetricType::Pointer _metaMetrics = MetaMetricType::New();
    OptimizerType::Pointer _optimizer = OptimizerType::New();

    QVector<ParametersType> params;
    for (unsigned int i = 0; i < _labelIndexes.size(); i++) {
        if (_labelIndexes[i].size() == 0) {
            continue;
        }
        cout << "# of indices: " << _labelIndexes[i].size() << endl;
        SubMetricType::Pointer metric = SubMetricType::New();
        InterpolatorType::Pointer interpolator = InterpolatorType::New();
        TransformType::Pointer transform = TransformType::New();
        metric->SetFixedImage(_fixedImage);
        metric->SetFixedImageIndexes(_labelIndexes[i]);
        metric->SetMovingImage(_movingImage);
        metric->SetInterpolator(interpolator);
        metric->SetTransform(dynamic_cast<SubMetricType::TransformType*>(transform.GetPointer()));
        metric->Initialize();
        _metaMetrics->AddMetric(metric);
        transform->SetFixedParameters(_centerOfRotation[i]);
        params.push_back(transform->GetParameters());
        cout << "Initial parameter: " << transform->GetParameters() << "\n";
        cout << "Center of rotation: " << transform->GetFixedParameters() << "\n";
    }

    MetaMetricType::ParametersType initialParams;
    initialParams.SetSize(_metaMetrics->GetNumberOfParameters());

    int k = 0;
    for (unsigned int i = 0; i < params.size(); i++) {
        for (int j = 0; j < params[i].GetSize(); j++) {
            initialParams[k++] = params[i][j];
        }
    }

    OptimizerType::ScalesType paramScales;
    paramScales.SetSize(_metaMetrics->GetNumberOfParameters());
    for (int i = 0; i < paramScales.GetSize(); i += 9) {
        paramScales[i] = 30;
        paramScales[i+1] = 30;
        paramScales[i+2] = 30;
        paramScales[i+3] = 0.5;
        paramScales[i+4] = 0.5;
        paramScales[i+5] = 0.5;
        paramScales[i+6] = 160;
        paramScales[i+7] = 160;
        paramScales[i+8] = 160;
    }
    _optimizer->SetScales(paramScales);
    _optimizer->SetCostFunction(_metaMetrics);
    _optimizer->SetInitialPosition(initialParams);
    _optimizer->AddObserver(IterationEvent(), this);
    _optimizer->StartOptimization();

    ParametersType param = _metaMetrics->GetParameters();
    _transformHistory.push_back(param);
}

LabelType::Pointer itkMultiRegMethod::TransformFixedLabel(int historyId) {
	return itk::SmartPointer<LabelType>(NULL);
}

void itkMultiRegMethod::WriteTransform(const char* filename, int historyId) {

}

