//
//  itkMyRegistration.h
//  itktools
//
//  Created by Joohwi Lee on 9/21/12.
//
//

#ifndef itktools_itkMyRegistration_h
#define itktools_itkMyRegistration_h

#include "vector"

#include "itkMyTypes.h"
#include "itkMyMetric.h"
#include "itkScaleVersor3DTransform.h"
#include "itkMyScaleVersor3DTransformOptimizer.h"
#include "itkCommand.h"
#include "itkTimer.h"
#include "itkRealTimeClock.h"
#include "QVector"

#define USE_GD_OPTIMIZER

template<class TTransform, class TOptimizer>
class RegistrationEngine : public itk::Command {
public:

    typedef RegistrationEngine Self;
    typedef itk::Command Superclass;
    typedef itk::SmartPointer<Self> Pointer;
    typedef itk::SmartPointer<const Self> ConstPointer;
    typedef itk::MyMetric<ImageType, ImageType> Metric;
    typedef TOptimizer OptimizerType;
    typedef QVector<typename TTransform::ParametersType> TransformHistoryType;


    itkNewMacro(RegistrationEngine);
    itkTypeMacro(RegistrationEngine, itk::Command);

private:
    ImageType::Pointer _movingImage;
    ImageType::Pointer _fixedImage;
    LabelType::Pointer _fixedImageLabel;

    ImageType::PointType _movingCenter;
    ImageType::PointType _fixedCenter;

    typename TTransform::Pointer _transform;
    typename TTransform::ParametersType _transformResult;
    typename TTransform::ParametersType _centerOfRotation;

    Metric::FixedImageIndexContainer _labelIndexes;

    char* _transformOut;
    char* _resampledOut;
    string _method;

    itk::RealTimeClock::Pointer _clock;
    itk::RealTimeClock::TimeStampType _lastTime;
    int _numOfIterations;
    int _updateInterval;
    bool _saveIntermediateImage;
    
    TransformHistoryType _transformHistory;
    TransformHistoryType _transformHistoryFixedParameters;
    
private:
    RegistrationEngine(const RegistrationEngine&);

protected:
    RegistrationEngine() {
        _updateInterval = 1;
        _clock = itk::RealTimeClock::New();
        _numOfIterations = 0;
    }
    
    virtual ~RegistrationEngine() {
    }


public:
    void Execute(itk::Object* caller, const itk::EventObject& event) {
        Execute((const Object*) caller, event);
    }

    void OnIterate(const OptimizerType* optimizer) {
        if (++ _numOfIterations % _updateInterval == 0) {
            itk::RealTimeClock::TimeStampType t = _clock->GetTimeInSeconds();
            std::cout << _numOfIterations << "\t" << optimizer->GetValue() << "\t" << optimizer->GetCurrentPosition() << "\t" << (t - _lastTime) << " secs" << std::endl;
            _transformHistory.push_back(optimizer->GetCurrentPosition());

            _lastTime = t;
        }
    }

    void Execute(const itk::Object* object, const itk::EventObject& event) {
        if (object == NULL) {
            std::cout << "Null sender is not processed..." << std::endl;
            return;
        }
        if (typeid(event) == typeid(itk::IterationEvent)) {
            const OptimizerType* optimizer = dynamic_cast<const OptimizerType*>(object);
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

    void Update() {
        this->Execute((const Object*) NULL, itk::IterationEvent());
    }
    

    TransformHistoryType GetTransformHistory() {
        return _transformHistory;
    }

    void WriteTransform(const char* filename, int historyId) {
        if (_transformHistory.size() == 0) {
            return;
        }

        while (historyId < 0) {
            historyId = historyId + _transformHistory.size();
        }
        
        while (historyId >= _transformHistory.size()) {
            historyId -= _transformHistory.size();
        }

        TransformWriterType::Pointer writer = TransformWriterType::New();
        writer->SetFileName(filename);
        typename TTransform::Pointer transform = TTransform::New();
        transform->SetParameters(_transformHistory.at(historyId));
        transform->SetFixedParameters(_transform->GetFixedParameters());
        writer->AddTransform(transform);
        writer->Update();
    }

    LabelType::Pointer TransformFixedLabel(int historyId) {
        if (_transformHistory.size() == 0 || historyId == 0) {
            return _fixedImageLabel;
        }

        typename TTransform::Pointer transform = TTransform::New();
        transform->SetParameters(_transformHistory.at(historyId - 1));
        transform->SetFixedParameters(_transform->GetFixedParameters());
        typename TTransform::InverseTransformBasePointer inverseTransform = transform->GetInverseTransform();
        typedef itk::ResampleImageFilter<LabelType, LabelType> ResampleFilter;
        ResampleFilter::Pointer resampler = ResampleFilter::New();
        try {
            if (_fixedImageLabel.IsNull()) {
                cout << "Input image is null ..." << endl;
                return _fixedImageLabel;
            }
            resampler->SetInput(_fixedImageLabel);
            resampler->SetTransform(inverseTransform);
            resampler->SetInterpolator(InterpolatorNN::New());
            resampler->SetReferenceImage(_fixedImageLabel);
            resampler->SetUseReferenceImage(true);
            resampler->Update();
        } catch (itk::ExceptionObject& ex) {
            cout << ex << endl;
        }
        return resampler->GetOutput();
    }

    void SetImages(ImageType::Pointer source, LabelType::Pointer sourceLabel, ImageType::Pointer target) {
        _movingImage = target;
        _fixedImage = source;
        _fixedImageLabel = sourceLabel;
        
        ImageType::SizeType szDst = _fixedImage->GetBufferedRegion().GetSize();
        itk::ContinuousIndex<double,3> szIdx;
        for (int i = 0; i < 3; i++) {
            szIdx[i] = szDst[i] / 2.0;
        }
        _fixedImage->TransformContinuousIndexToPhysicalPoint(szIdx, _fixedCenter);
    }

 
    void RunRegistration() {
        _transform = TTransform::New();
        Metric::Pointer metric = Metric::New();
        metric->SetFixedImage(_fixedImage);
        bool useIndexes = false;
        if (useIndexes) {
            _centerOfRotation.SetSize(ImageType::ImageDimension);
            for (int i = 0; i < ImageType::ImageDimension; i++) {
                _centerOfRotation[i] = i;
            }
            itk::ImageRegionConstIteratorWithIndex<LabelType> labelIter(_fixedImageLabel, _fixedImageLabel->GetBufferedRegion());
            int nPixels = 0;
            for (labelIter.GoToBegin(); !labelIter.IsAtEnd(); ++labelIter) {
                LabelType::PixelType label = labelIter.Get();
                if (label > 0) {
                    _labelIndexes.push_back(labelIter.GetIndex());
                    for (int i = 0; i < ImageType::ImageDimension; i++) {
                        _centerOfRotation[i] += labelIter.GetIndex()[i];
                    }
                    nPixels ++;
                }
            }
            for (int i = 0; i < ImageType::ImageDimension; i++) {
                _centerOfRotation[i] /= nPixels;
            }
            metric->SetFixedImageIndexes(_labelIndexes);
            _transform->SetFixedParameters(_centerOfRotation);
        } else {
            metric->SetFixedImageRegion(_fixedImage->GetBufferedRegion());
            metric->SetUseAllPixels(true);
            _centerOfRotation.SetSize(ImageType::ImageDimension);
            for (int i = 0; i < 3; i++) {
                _centerOfRotation[i] = _fixedCenter[i];
            }
            _transform->SetFixedParameters(_centerOfRotation);
        }

        metric->SetMovingImage(_movingImage);
        metric->SetInterpolator(NNImageInterpolator::New());
        metric->SetTransform(_transform);
        metric->Initialize();

        typename OptimizerType::Pointer opti = OptimizerType::New();
        opti->SetCostFunction(metric);
        typename OptimizerType::ScalesType scales;
        scales.SetSize(TTransform::ParametersDimension);
        scales.Fill(1);
        scales[0] = scales[1] = scales[2] = 30;
        scales[3] = scales[4] = scales[5] = 0.05;
        scales[6] = scales[7] = scales[8] = 160;

        opti->SetScales(scales);

        const int maxIters = 50;
#ifdef USE_CG_Optimizer
        opti->SetMaximumIteration(maxIters);
        opti->SetMaximumLineIteration(10);
        opti->SetUseUnitLengthGradient(true);
        opti->SetStepLength(1);
        opti->SetToFletchReeves();
#endif

#ifdef USE_GD_OPTIMIZER
        opti->SetNumberOfIterations(maxIters);
        opti->SetMinimumStepLength(1e-4);
        opti->SetMaximumStepLength(1);
        opti->SetRelaxationFactor(.5);
        opti->SetGradientMagnitudeTolerance(1e-4);
#endif

        opti->SetInitialPosition(_transform->GetParameters());
        opti->AddObserver(itk::StartEvent(), this);
        opti->AddObserver(itk::IterationEvent(), this);
        opti->StartOptimization();
        cout << "Current Cost: " << opti->GetValue() << endl;
        _transformResult = opti->GetCurrentPosition();
        _transform->SetParameters(opti->GetCurrentPosition());
    }

    void SaveTransform() {
        TransformWriterType::Pointer writer = TransformWriterType::New();
        //writer->AddTransform(_initialTransform);
        writer->AddTransform(_transform);
        writer->SetFileName(_transformOut);
        writer->Update();
    }


};


typedef RegistrationEngine<itk::ScaleVersor3DTransform<double>, itk::MyScaleVersor3DTransformOptimizer > ScaleRegistration;


#endif
