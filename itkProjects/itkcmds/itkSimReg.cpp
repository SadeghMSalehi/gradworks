/*
 * itkSimReg.cpp
 *
 *  Created on: Sep 5, 2012
 *      Author: joohwi
 */

#include "itkImageIO.h"
#include "itkExceptionObject.h"
#include "itkAffineTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkImageRegistrationMethod.h"
#include "itkResampleImageFilter.h"
#include "itkContinuousIndex.h"
#include "itkCropImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkLabelGeometryImageFilter.h"
#include "itkMath.h"
#include "itkTimer.h"
#include "iostream"
#include "fstream"
#include "vector"
#include "itkMyMetric.h"
#include "itkMyFRPROptimizer.h"
#include "MyFunction.h"
#include "MatrixCode.h"
#include "itkCommand.h"
#include "itkSingleValuedNonLinearOptimizer.h"
#include "itkRealTimeClock.h"
#include "itksys/hash_map.hxx"
#include "itkOptimizationReporter.h"
#include "itkExtractImageFilter.h"
#include "itkVersorRigid3DTransform.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkSimilarity3DTransform.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkMathCode.h"
#include "itkFastMeanSquaredMetric.h"
#include "itkMeanSquaresImageToImageMetric.h"

class SimilarityRegistration : public itk::Command {
public:
    typedef itk::Image<double,3> ImageType;
    typedef itk::Image<unsigned short, 3> LabelType;
    typedef itk::LabelGeometryImageFilter<LabelType> LabelFilterType;
    typedef itk::ExtractImageFilter<ImageType, ImageType> CropFilterType;
    typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;
    typedef itk::NearestNeighborInterpolateImageFunction<ImageType> NearestNeighborInterpolatorType;
    typedef itk::Similarity3DTransform<double> TransformType;
    typedef itk::ImageRegistrationMethod<ImageType, ImageType> RegistrationType;
    typedef itk::MyFRPROptimizer OptimizerType;
    typedef itk::MeanSquaresImageToImageMetric<ImageType, ImageType> MetricType;
    typedef itk::ExceptionObject itkException;
    typedef itk::ResampleImageFilter<ImageType, ImageType> ResamplerType;

    typedef SimilarityRegistration Self;
    typedef itk::Command Superclass;
    typedef itk::SmartPointer<Self> Pointer;

    itkTypeMacro(SimilarityRegistration, itk::Command);
    itkNewMacro(SimilarityRegistration);

private:
    ImageType::Pointer _srcImg;
    ImageType::Pointer _srcROIImg;
    ImageType::Pointer _dstImg;
    LabelType::Pointer _labelImg;
    MetricType::Pointer _metric;
    OptimizerType::Pointer _optimizer;
    TransformType::Pointer _transform;
    RegistrationType::TransformPointer _finalTransform;
    InterpolatorType::Pointer _interpolator;

    int _iteration;
    int _updateInterval;

public:
    void OnStartEvent(const itk::Object* caller) {
        cout << "Registration start ..." << endl;
        _iteration = 0;
        _updateInterval = 1;
    }

    void OnEndEvent(const itk::Object* caller) {
        cout << "Registration end ..." << endl;
    }

    void OnIterationEvent(const itk::Object* caller) {
        if (_iteration % _updateInterval != 0) {
            return;
        }
        cout << _optimizer->GetCurrentCost() << "\t" << _transform->GetParameters() << "\t" << _transform->GetFixedParameters() << endl;

        char buffName[128];
        sprintf(buffName, "resampled_%02d.nrrd", _optimizer->GetCurrentIteration());
        ResampleResult(buffName);
    }

    void Execute(itk::Object* caller, const itk::EventObject& event) {
        Execute((const itk::Object*) caller, event);
    }

    void Execute(const itk::Object* caller, const itk::EventObject& event) {
        if (typeid(itk::StartEvent) == typeid(event)) {
            OnStartEvent(caller);
        } else if (typeid(itk::EndEvent) == typeid(event)) {
            OnEndEvent(caller);
        } else if (typeid(itk::IterationEvent) == typeid(event)) {
            OnIterationEvent(caller);
        }
    }


    void Update() {
        this->Execute((const itk::Object*) NULL, itk::IterationEvent());
    }


    void LoadFiles(const char* name, const char* labelname, const char* dstname) {
        itkcmds::itkImageIO<ImageType> imageIO;
        itkcmds::itkImageIO<LabelType> labelIO;

        _srcImg = imageIO.ReadImageT(name);
        _labelImg = labelIO.ReadImageT(labelname);
        _dstImg = imageIO.ReadImageT(dstname);
    }

    void SetupMetrics() {
        _metric = MetricType::New();
        _transform = TransformType::New();
        _metric->SetFixedImage(_srcImg);
        _metric->SetMovingImage(_dstImg);

        MetricType::FixedImageIndexContainer indexes;
        itkcmds::ComputeLabelIndexes<LabelType, MetricType>(_labelImg, indexes);

        cout << "# of indexes: " << indexes.size() << endl;
        _metric->SetFixedImageIndexes(indexes);

        _metric->SetInterpolator(InterpolatorType::New());
        _metric->SetTransform(_transform);
        _metric->Initialize();
    }

    int RunRegistration() {
        _optimizer = OptimizerType::New();
        
        _optimizer->AddObserver(StartEvent(), this);
        _optimizer->AddObserver(IterationEvent(), this);
        _optimizer->AddObserver(EndEvent(), this);

        OptimizerType::ScalesType scales;
        scales.SetSize(_transform->GetNumberOfParameters());
        scales[0] = scales[1] = scales[2] = 10;
        scales[3] = scales[4] = scales[5] = 0.1;
        scales[6] = 1;
        //_optimizer->SetScales(scales);
        //_optimizer->SetToFletchReeves();
        _optimizer->SetCostFunction(_metric);
        TransformType::ParametersType initialParameter = _transform->GetParameters();
        _optimizer->SetInitialPosition(initialParameter);
        cout << "Center of rotation: " << _transform->GetFixedParameters() << endl;
        _optimizer->StartOptimization();
        cout << _optimizer->GetCurrentPosition() << endl;
        return EXIT_SUCCESS;
    }

    void ResampleResult(const char* name) {
        NearestNeighborInterpolatorType::Pointer resampleInterpolator = NearestNeighborInterpolatorType::New();
        ResamplerType::Pointer resample = ResamplerType::New();
        resample->SetTransform(_transform);
        resample->SetInput(_srcImg);
        resample->UseReferenceImageOn();
        resample->SetReferenceImage(_dstImg);
        resample->SetInterpolator(resampleInterpolator);
        resample->Update();
        ImageType::Pointer resampledImage = resample->GetOutput();
        itkcmds::itkImageIO<ImageType> io;
        io.WriteImageT(name, resampledImage);
    }

    void Run(int argc, char* argv[]) {
        //        TransformType::Pointer transform = TransformType::New();
        //        TransformType::ParametersType params;
        //        params.SetSize(transform->GetNumberOfParameters());
        //        params.Fill(0);
        //        params[0] = 0;
        //        params[1] = 0;
        //        params[2] = 1*sin(itk::Math::pi_over_4);
        //        params[3] = 0;
        //        params[4] = 0;
        //        params[5] = 0;
        //
        //        cout << "Parameters: " << params << endl;
        //        transform->SetParameters(params);
        //        TransformType::InputPointType p;
        //        p[0] = 1;
        //        p[1] = 0;
        //        p[2] = 0;
        //        TransformType::OutputPointType q = transform->TransformPoint(p);
        //        cout << p << "=>" << q << endl;

        LoadFiles(argv[1], argv[2], argv[3]);
        SetupMetrics();
        RunRegistration();
        //ResampleResult(argv[4]);
    }

    void Test(int argc, char* argv[]) {
        LoadFiles(argv[1], argv[2], argv[3]);

        typedef FastMeanSquaredMetric<ImageType, LabelType, TransformType, InterpolatorType> FastMetricType;
        FastMetricType::Pointer metric = FastMetricType::New();

        metric->SetMovingImage(_srcImg);
        metric->SetMovingImageMask(_labelImg);
        metric->SetFixedImage(_dstImg);

        TransformType::Pointer transform = TransformType::New();
        metric->SetTransform(transform);
        metric->Initialize();

        TransformType::ParametersType initialParam = transform->GetParameters();
        initialParam[2] = 1*sin(itk::Math::pi_over_4);
        initialParam[3] = 0;
        initialParam[4] = 0;
        initialParam[5] = 0;
        initialParam[6] = 1;

        TransformType::ParametersType fixedParam = transform->GetFixedParameters();
        fixedParam[0] = 0;
        fixedParam[1] = 0;
        fixedParam[2] = 0;
        transform->SetFixedParameters(fixedParam);
        transform->SetParameters(initialParam);

        typedef itkMathCode<ImageType, TransformType> MathCodeType;
        MathCodeType mathCode;
        MathCodeType::Matrix imageTransform;
        mathCode.createIndex2IndexTransformMatrix(_srcImg, transform, _dstImg, imageTransform);
        //mathCode.convertTransformToMatrix(transform, imageTransform);

        TransformType::InputPointType point;
        point[0] = 3;
        point[1] = 7;
        point[2] = 9;
        TransformType::OutputPointType qPoint = transform->TransformPoint(point);

        MathCodeType::Vector x(point[0],point[1],point[2],1), y;
        MathCode::mult(imageTransform, x, y);

        cout << "Transformed1: " << point << " => " << qPoint << endl;
        cout << "Transformed2: " << x << " => " << y << endl;
        //metric->GetValue(initialParam);
        
    }
private:
    
};

int main(int argc, char* argv[]) {
    SimilarityRegistration::Pointer engine = SimilarityRegistration::New();
    //engine->Test(argc, argv);
    engine->Run(argc, argv);
}
