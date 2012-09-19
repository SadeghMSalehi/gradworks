#include "itkImageIO.h"
#include "itkTimer.h"
#include "itkExceptionObject.h"
#include "itkSimilarity2DTransform.h"
#include "itkSimilarity3DTransform.h"
#include "itkAffineTransform.h"
#include "itkLinearInterpolateImageFunction.h"

#include "itkImageRegistrationMethod.h"
#include "itkResampleImageFilter.h"
#include "itkContinuousIndex.h"
#include "itkMath.h"
#include "iostream"
#include "itkMyFRPROptimizer.h"
#include "itkMetaMetrics.h"
#include "itkOptimizationReporter.h"
#include "itkMyMetric.h"

#include "itkCropImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkResampler.h"

using namespace std;
using namespace itk;

typedef itk::ExceptionObject itkException;
typedef itk::Image<double, 3> ImageType;
typedef itk::Image<unsigned short, 3> LabelType;
typedef itk::Similarity3DTransform<double> TransformType;
//typedef itk::AffineTransform<double, 3> TransformType;
typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;

typedef itk::MyFRPROptimizer OptimizerType;
typedef itk::MyMetric<ImageType, ImageType> SubMetricType;
typedef MetaMetrics<SubMetricType> MetaMetricType;

#define angle2rad(a) (a*itk::Math::pi/180.0)

class MultipleRegionRegistration {
private:
    typedef vector<ImageType::Pointer> ImageVector;
    MetaMetricType::Pointer _metaMetrics;
    MetaMetricType::ParametersType _currentParams;
    ImageType::Pointer _fixedImage;
    ImageType::Pointer _movingImage;
    LabelType::Pointer _fixedLabel;
    itkcmds::itkImageIO<ImageType> _imageIO;
    itkcmds::itkImageIO<LabelType> _labelIO;
    OptimizerType::Pointer _optimizer;
    OptimizationReporter<OptimizerType>::Pointer _reporter;
    vector< SubMetricType::FixedImageIndexContainer> _labelIndices;
    typedef itk::ResampleImageFilter<ImageType, ImageType> ResamplerType;

public:
    void LoadImages(int argc, char* argv[]) {
        _fixedImage = _imageIO.ReadImageT(argv[1]);
        _fixedLabel = _labelIO.ReadImageT(argv[2]);
        _movingImage = _imageIO.ReadImageT(argv[3]);
    }

    void Initialize() {
        _metaMetrics = MetaMetricType::New();
        _labelIndices.resize(50);

        Timer timer;
        timer.start();
        ImageRegionConstIteratorWithIndex<LabelType> labelIter(_fixedLabel, _fixedLabel->GetBufferedRegion());
        for (labelIter.GoToBegin(); !labelIter.IsAtEnd(); ++labelIter) {
            LabelType::PixelType label = labelIter.Get();
            if (label > 0) {
                _labelIndices[label].push_back(labelIter.GetIndex());
            }
        }
        timer.stop();
        cout << "Label Indexing Time: " << timer.getElapsedTimeInMilliSec() << " ms" << endl;

        _reporter = OptimizationReporter<OptimizerType>::New();
        _optimizer = OptimizerType::New();

        vector<TransformType::ParametersType> params;
        for (unsigned int i = 0; i < _labelIndices.size(); i++) {
            if (_labelIndices[i].size() == 0) {
                continue;
            }
            cout << "# of indices: " << _labelIndices[i].size() << endl;
            SubMetricType::Pointer metric = SubMetricType::New();
            InterpolatorType::Pointer interpolator = InterpolatorType::New();
            TransformType::Pointer transform = TransformType::New();
            metric->SetFixedImage(_fixedImage);
            metric->SetFixedImageIndexes(_labelIndices[i]);
            metric->SetMovingImage(_movingImage);
            metric->SetInterpolator(interpolator);
            metric->SetTransform(transform);
            metric->Initialize();
            _metaMetrics->AddMetric(metric);
            cout << "Initial parameter: " << transform->GetParameters() << endl;
            cout << "Center of rotation: " << transform->GetFixedParameters() << endl;
            params.push_back(transform->GetParameters());
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
        for (int i = 0; i < paramScales.GetSize(); i += 7) {
            paramScales[i+6] = 100;
            paramScales[i] = 1;
            paramScales[i+1] = 1;
            paramScales[i+2] = 1;
            paramScales[i+3] = 1;
            paramScales[i+4] = 1;
            paramScales[i+5] = 1;
        }
        _optimizer->SetMaximumIteration(100);
        _optimizer->SetScales(paramScales);
        _optimizer->SetCostFunction(_metaMetrics);
        _optimizer->SetInitialPosition(initialParams);
        _optimizer->SetUseUnitLengthGradient(true);
        _optimizer->SetStepLength(0.25);
        _optimizer->AddObserver(StartEvent(), _reporter);
        _optimizer->AddObserver(IterationEvent(), _reporter);
        _optimizer->AddObserver(EndEvent(), _reporter);
    }

    void StartOptimization() {
        _optimizer->StartOptimization();
    }

    MetaMetricType::ParametersType GetCurrentParameter() {
        return _optimizer->GetCurrentPosition();
    }

    void ResampleOutput(int argc, char* argv[]) {
        itkResampler<ImageType, LabelType, TransformType> resampler;

        TransformType::Pointer transform = TransformType::New();
        transform->SetParameters(_optimizer->GetCurrentPosition());
        transform->SetFixedParameters(_metaMetrics->GetMetricPointer(0)->GetFixedParameters());
        resampler.SetSourceMask(_fixedLabel);

        resampler.SetTransform(transform);
        resampler.ResampleBackward();
        _imageIO.WriteImageT(argv[4], resampler.GetOutput());
    }
};

int main(int argc, char* argv[]) {
    MultipleRegionRegistration regEngine;
    regEngine.LoadImages(argc, argv);
    regEngine.Initialize();
    regEngine.StartOptimization();
    return 0;
}