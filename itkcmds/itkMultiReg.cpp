#include "itkImageCommon.h"
#include "itkExceptionObject.h"
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

using namespace std;
using namespace itk;

typedef itk::ExceptionObject itkException;
typedef itk::Image<double, 3> ImageType;
typedef itk::Image<short, 3> LabelType;
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
    ImageVector _movingImages;
    itkcmds::itkImageIO<ImageType> _imageIO;
    OptimizerType::Pointer _optimizer;
    OptimizationReporter::Pointer _reporter;

public:
    void loadImages(int argc, char* argv[]) {
        ImageType::Pointer fixed = _imageIO.ReadImageT(argv[1]);
        for (int i = 2; i < argc; i++) {
            ImageType::Pointer moving = _imageIO.ReadImageT(argv[i]);
            _movingImages.push_back(moving);
        }
        _reporter = OptimizationReporter::New();
    }

    void initialize() {
        _metaMetrics = MetaMetricType::New();
        SubMetricType::ParametersType initialParams;
        initialParams.SetSize(_metaMetrics->GetNumberOfParameters());
        int nOffset = 0;
        for (ImageVector::size_type i = 0; i < _movingImages.size(); i++) {
            SubMetricType::Pointer metric = SubMetricType::New();
            InterpolatorType::Pointer interpolator = InterpolatorType::New();
            TransformType::Pointer transform = TransformType::New();
            metric->SetFixedImage(_fixedImage);
            metric->SetMovingImage(_movingImages[i]);
            metric->SetFixedImageRegion(_fixedImage->GetBufferedRegion());
            metric->UseAllPixelsOn();
            metric->SetInterpolator(interpolator);
            metric->SetTransform(dynamic_cast<SubMetricType::TransformType*>(transform.GetPointer()));
            metric->Initialize();
            _metaMetrics->AddMetric(metric);
            int nParam = transform->GetNumberOfParameters();
            TransformType::ParametersType subParam = transform->GetParameters();
            for (int i = 0; i < nParam; i++) {
                initialParams[nOffset + i] = subParam[i];
            }
        }
        _optimizer->SetCostFunction(_metaMetrics);
        _optimizer->SetInitialPosition(initialParams);
        _optimizer->AddObserver(StartEvent(), _reporter);
        _optimizer->AddObserver(IterationEvent(), _reporter);
        _optimizer->AddObserver(EndEvent(), _reporter);
    }

    void startOptimization() {
        _optimizer->StartOptimization();
    }

    MetaMetricType::ParametersType getCurrentParameter() {
        _optimizer->GetCurrentPosition();
    }
};

int mainRegistration(int argc, char* argv[]) {
    MultipleRegionRegistration regEngine;
    regEngine.loadImages(argc, argv);
    regEngine.initialize();
    regEngine.startOptimization();
    return 0;
}

void computeMSE(ImageType::Pointer src, ImageType::Pointer dst) {
	double* srcBuff = src->GetBufferPointer();
	double* dstBuff = dst->GetBufferPointer();

	ImageType::RegionType srcRegion = src->GetBufferedRegion();
	int nPixels = srcRegion.GetNumberOfPixels();

	double mse = 0;
	int nCount = 0;
	for (int i = 0; i < nPixels; i++) {
		double diff = srcBuff[i] - dstBuff[i];
		mse += diff * diff;
		nCount++;
	}
	mse /= nPixels;
	cout << "Mean Squared Error: " << mse << "; # pixels: " << nCount << endl;
}

int mainCrop(int argc, char* argv[]) {
    typedef itk::CropImageFilter<ImageType,ImageType> CropFilterType;
    SubMetricType::Pointer metric = SubMetricType::New();
    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    TransformType::Pointer transform = TransformType::New();

    if (argc < 3) {
        cout << "multiReg fixedImage movingImage" << endl;
        return 1;
    }
    itkcmds::itkImageIO<ImageType> imageIO;
    ImageType::Pointer fixed = imageIO.ReadImageT(argv[1]);
    ImageType::Pointer moving = imageIO.ReadImageT(argv[2]);

    metric->SetInterpolator(interpolator);
    metric->SetTransform(dynamic_cast<SubMetricType::TransformType*> (transform.GetPointer()));
    metric->UseAllPixelsOn();
    metric->SetFixedImage(fixed);
    metric->SetFixedImageRegion(fixed->GetBufferedRegion());
    metric->SetMovingImage(moving);
    metric->Initialize();

    transform->SetIdentity();
    
    SubMetricType::ParametersType initialParam = transform->GetParameters();
    cout << "Initial Parameter: " << initialParam << endl;
    computeMSE(fixed, moving);
    cout << "Initial value: " << metric->GetValue(initialParam) << endl;

    OptimizationReporter::Pointer reporter = OptimizationReporter::New();

    OptimizerType::Pointer optimizer = OptimizerType::New();
    optimizer->SetCostFunction(metric);
    optimizer->SetInitialPosition(initialParam);
    optimizer->AddObserver(itk::StartEvent(), reporter);
    optimizer->AddObserver(itk::IterationEvent(), reporter);
    optimizer->AddObserver(itk::EndEvent(), reporter);



    try {
    optimizer->StartOptimization();
    } catch (itk::ExceptionObject& ex) {
        cout << "Exception caught: " << ex << endl;
    }
    TransformType::ParametersType solution = optimizer->GetCurrentPosition();
    cout << solution << endl;
    return 0;
}

int main(int argc, char* argv[]) {
    return mainCrop(argc, argv);
}
