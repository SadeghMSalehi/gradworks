/*
 * itkMultiRegMethod.h
 *
 *  Created on: Sep 24, 2012
 *      Author: joohwi
 */

#ifndef ITKMULTIREGMETHOD_H_
#define ITKMULTIREGMETHOD_H_

#include "itkCommand.h"
#include "itkMyRegularStepGradientDescentOptimizer.h"
#include "itkMyScaleVersor3DTransformOptimizer.h"
#include "itkMyTypes.h"
#include "itkMetaMetrics.h"
#include "itkMyMetric.h"
#include "itkRealTimeClock.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkScaleVersor3DTransform.h"
#include "QVector"



class itkMultiRegMethod: public itk::Command {
public:
	typedef itkMultiRegMethod Self;
	typedef itk::Command Superclass;
	typedef itk::SmartPointer<Self> Pointer;
	typedef itk::SmartPointer<const Self> ConstPointer;
	typedef itk::MyMetric<ImageType, ImageType> Metric;
	typedef itk::LinearInterpolateImageFunction<ImageType> InterpolatorType;
	typedef itk::MyRegularStepGradientDescentOptimizer OptimizerType;
	typedef itk::ScaleVersor3DTransform<double> TransformType;
	typedef itk::SingleValuedCostFunction::ParametersType ParametersType;
	typedef QVector<ParametersType> TransformHistoryType;

	itkNewMacro(itkMultiRegMethod)
	;itkTypeMacro(itkMultiRegMethod, itk::Command)
	;

private:
	ImageType::Pointer _movingImage;
	ImageType::Pointer _fixedImage;
	LabelType::Pointer _fixedImageLabel;

	ImageType::PointType _movingCenter;

	QVector<ImageType::PointType> _fixedCenter;
	QVector<typename TransformType::Pointer> _transform;
	QVector<ParametersType> _transformResult;
	QVector<ParametersType> _centerOfRotation;
	QVector<Metric::FixedImageIndexContainer> _labelIndexes;

	std::string _method;

	itk::RealTimeClock::Pointer _clock;
	int _numOfIterations;
	int _updateInterval;
	itk::RealTimeClock::TimeStampType _lastTime;

	TransformHistoryType _transformHistory;
	TransformHistoryType _transformHistoryFixedParameters;

private:
	itkMultiRegMethod(const itkMultiRegMethod&);

protected:
	itkMultiRegMethod() {
		_updateInterval = 1;
		_clock = itk::RealTimeClock::New();
		_numOfIterations = 0;
		_lastTime = 0;
	}

	virtual ~itkMultiRegMethod() {
	}

public:
	void Execute(itk::Object* caller, const itk::EventObject& event);
	void OnIterate(const OptimizerType* optimizer);
	void Execute(const itk::Object* object, const itk::EventObject& event);
	void Update();

public:

	void SetImages(ImageType::Pointer source, LabelType::Pointer sourceLabel,
			ImageType::Pointer target);
	void RunRegistration();
	LabelType::Pointer TransformFixedLabel(int historyId);
	void WriteTransform(const char* filename, int historyId);

	TransformHistoryType GetTransformHistory() {
		return _transformHistory;
	}

};


typedef itkMultiRegMethod RegistrationMethod;

#endif /* ITKMULTIREGMETHOD_H_ */
