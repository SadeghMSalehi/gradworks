/*
 * itkImageReg.cpp
 *
 *  My own registration program with 9 parameters (3*scale, 3*rotation, 3*translation)
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
#include "itkScaleVersor3DTransform.h"

using namespace std;
using namespace itk;

#define angle2rad(a) (a*itk::Math::pi/180.0)

class RegionalRegistration : public itk::Command {
public:
    typedef itk::Image<double,3> ImageType;
    typedef itk::Image<unsigned short, 3> LabelType;
    typedef itk::LabelGeometryImageFilter<LabelType> LabelFilterType;
    typedef itk::ExtractImageFilter<ImageType, ImageType> CropFilterType;
    typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;
    typedef itk::NearestNeighborInterpolateImageFunction<ImageType> NearestNeighborInterpolatorType;

    typedef itk::Rigid3DTransform<double> Transform1Type;
    typedef itk::ScaleVersor3DTransform<double> Transform2Type;

    typedef itk::ImageRegistrationMethod<ImageType, ImageType> RegistrationType;
    typedef itk::MyFRPROptimizer OptimizerType;
    typedef itk::MyMetric<ImageType, ImageType> MetricType;
    typedef itk::ExceptionObject itkException;
    typedef itk::ResampleImageFilter<ImageType, ImageType> ResamplerType;

    typedef RegionalRegistration Self;
    typedef itk::Command Superclass;
    typedef itk::SmartPointer<Self> Pointer;

    itkTypeMacro(RegionalRegistration, itk::Command);
    itkNewMacro(RegionalRegistration);

private:
    ImageType::Pointer _srcImg;
    ImageType::Pointer _srcROIImg;
    ImageType::Pointer _dstImg;
    LabelType::Pointer _labelImg;
    MetricType::Pointer _metric;
    OptimizerType::Pointer _optimizer;
    Transform1Type::Pointer _transform1;
    Transform2Type::Pointer _transform2;
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
//        cout << _optimizer->GetCurrentCost() << "\t" << _transform->GetParameters() << "\t" << _transform->GetFixedParameters() << endl;
//
//        char buffName[128];
//        sprintf(buffName, "resampled_%02d.nrrd", _optimizer->GetCurrentIteration());
//        ResampleResult(buffName);
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
        LabelFilterType::Pointer labelFilter = LabelFilterType::New();
        labelFilter->SetInput(_labelImg);
        labelFilter->CalculateOrientedBoundingBoxOff();
        labelFilter->Update();

        LabelFilterType::BoundingBoxType bboxCoords = labelFilter->GetBoundingBox(3);
        ImageType::RegionType cropRegion;
        ImageType::IndexType lowerIndex, upperIndex;
        for (int i = 0; i < 6; i+=2) {
            lowerIndex[i/2] = bboxCoords[i];
            upperIndex[i/2] = bboxCoords[i+1];
        }
        cropRegion.SetIndex(lowerIndex);
        cropRegion.SetUpperIndex(upperIndex);

        CropFilterType::Pointer cropFilter = CropFilterType::New();
        cropFilter->SetInput(_srcImg);
        cropFilter->SetExtractionRegion(cropRegion);
        cropFilter->Update();
        _srcImg = cropFilter->GetOutput();

        itkcmds::itkImageIO<ImageType> imageIO;
        imageIO.WriteImageT("crop_src.nrrd", _srcImg);
        
        _metric = MetricType::New();
        _transform1 = Transform1Type::New();
        _interpolator = InterpolatorType::New();
        _optimizer = OptimizerType::New();
        _metric->SetUseAllPixels(true);
        _metric->SetFixedImage(_dstImg);
        _metric->SetFixedImageRegion(_dstImg->GetBufferedRegion());
        _metric->SetMovingImage(_srcImg);
        _metric->SetInterpolator(_interpolator);
        _metric->SetTransform(_transform1);
    }

    int RunRegistration() {
        _optimizer->AddObserver(StartEvent(), this);
        _optimizer->AddObserver(IterationEvent(), this);
        _optimizer->AddObserver(EndEvent(), this);

        OptimizerType::ScalesType scales;
        scales.SetSize(_transform->GetNumberOfParameters());
        scales[0] = scales[1] = scales[2] = 10;
        scales[3] = scales[4] = scales[5] = 0.1;
        scales[6] = 1;
        _optimizer->SetScales(scales);
        _optimizer->SetToFletchReeves();

        RegistrationType::Pointer registration = RegistrationType::New();
        registration->SetMetric(dynamic_cast<RegistrationType::MetricType*>(_metric.GetPointer()));
        registration->SetOptimizer(_optimizer);
        registration->SetTransform(
                                   static_cast<RegistrationType::TransformPointer>(_transform));
        registration->SetInterpolator(_interpolator);
        registration->SetFixedImage(_dstImg);
        registration->SetMovingImage(_srcImg);
        registration->SetFixedImageRegion(_srcImg->GetBufferedRegion());
        TransformType::ParametersType initialParameter = _transform->GetParameters();
        
        cout << "Initial parameter: " << initialParameter << endl;
        registration->SetInitialTransformParameters(initialParameter);
        
        try {
            registration->Update();
        }
        catch (itkException &ex) {
            std::cerr << ex;
            std::cerr << std::endl;
            return EXIT_FAILURE;
        }

        _finalTransform = registration->GetTransform();
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
        ResampleResult(argv[4]);
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
    RegionalRegistration::Pointer engine = RegionalRegistration::New();
    //engine->Test(argc, argv);
    engine.Run(argc, argv);
}


/*

 typedef itk::ExceptionObject itkException;
 typedef itk::Image<double, 3> ImageType;
 typedef itk::Image<unsigned short, 3> LabelType;
 typedef itk::AffineTransform<double, 3> TransformType;
 typedef itk::MyFRPROptimizer OptimizerType;
 typedef itk::MyMetric<ImageType, ImageType> MetricType;
 typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;
 typedef itk::ImageRegistrationMethod<ImageType, ImageType> RegistrationType;
 typedef itk::ResampleImageFilter<ImageType, ImageType> ResamplerType;
 typedef itk::CropImageFilter<ImageType,ImageType> CropFilterType;
 typedef itk::LabelGeometryImageFilter<ImageType,ImageType> LabelGeometryFilter;
void createMat4FromMat(const ImageType::DirectionType &mat,
                       MathCode::Mat4& matOut) {
	for (int r = 0; r < 3; r++) {
		for (int c = 0; c < 3; c++) {
			matOut[r][c] = mat[r][c];
		}
	}
}

void createMat4FromVec(const ImageType::SpacingType &spacing,
                       MathCode::Mat4& matOut, bool inverse =  false) {
	for (int i = 0; i < 3; i++) {
		matOut[i][i] = spacing[i];
		if (inverse) {
			matOut[i][i] = 1 / matOut[i][i];
		}
	}
}

void createIndex2IndexTransformMatrix(ImageType::Pointer src,
                                      ImageType::Pointer ref, TransformType::MatrixType &transformMatrix,
                                      MathCode::Mat4& matOut) {
	MathCode::Mat4 srcDirection, srcSpacing, srcImageToWorld;
	createMat4FromMat(src->GetDirection(), srcDirection);
	createMat4FromVec(src->GetSpacing(), srcSpacing);
	MathCode::mult(srcDirection, srcSpacing, srcImageToWorld);
	for (int i = 0; i < 3; i++) {
		srcImageToWorld[i][3] = src->GetOrigin()[i];
	}

	MathCode::Mat4 srcTransform, srcToRefWorld;
	createMat4FromMat(static_cast<ImageType::DirectionType>(transformMatrix),
                      srcTransform);
	mult(srcTransform, srcImageToWorld, srcToRefWorld);

	MathCode::Mat4 refDirection, refSpacing, refSpacingDirection,
    refWorldToImage;
	createMat4FromVec(ref->GetSpacing(), refSpacing, true);
	createMat4FromMat(ref->GetInverseDirection(), refDirection);

	mult(refSpacing, refDirection, refSpacingDirection);
	for (int i = 0; i < 3; i++) {
		srcToRefWorld[i][3] -= ref->GetOrigin()[i];
	}
	cout << "Src To Ref World: " << srcToRefWorld << endl;
	cout << "Ref Spacing: " << refSpacing << endl;
	cout << "Ref SpacingDirection: " << refSpacingDirection << endl;
	mult(refSpacingDirection, srcToRefWorld, refWorldToImage);
	matOut.copyFrom(refWorldToImage);
}

ImageType::Pointer resampleImageRaw(ImageType::Pointer src,
                                    ImageType::Pointer ref, TransformType::Pointer transform) {
	MathCode::Mat4 srcIndexToRefIndex, refIndexToSrcIndex;
	TransformType::MatrixType transformMatrix = transform->GetMatrix();
	createIndex2IndexTransformMatrix(src, ref, transformMatrix,
                                     srcIndexToRefIndex);

	ImageType::RegionType srcRegion = src->GetBufferedRegion();
	double* srcBuff = src->GetBufferPointer();
	ImageType::RegionType refRegion = ref->GetBufferedRegion();
	ImageType::SizeType refSize = refRegion.GetSize();

	itkcmds::itkImageIO<ImageType> imageIO;
	ImageType::Pointer dst = imageIO.NewImageT(src);
	dst->FillBuffer(0);
	double* dstBuff = dst->GetBufferPointer();
	srcIndexToRefIndex.inverse(refIndexToSrcIndex);

	InterpolatorType::Pointer interp = InterpolatorType::New();
	interp->SetInputImage(src);

	for (int z = 0; z < refSize[2]; z++) {
		for (int y = 0; y < refSize[1]; y++) {
			for (int x = 0; x < refSize[0]; x++) {
				MathCode::Vec4 dstPoint(x, y, z, 1);
				MathCode::Vec4 srcPoint;
				mult(refIndexToSrcIndex, dstPoint, srcPoint);
				itk::ContinuousIndex<double,3> srcIndex;
				srcIndex[0] = srcPoint[0];
				srcIndex[1] = srcPoint[1];
				srcIndex[2] = srcPoint[2];
				if (!srcRegion.IsInside(srcIndex)) {
					continue;
				}

			}
		}
	}

	return src;
}

ImageType::Pointer resampleImageITK(ImageType::Pointer src,
                                    TransformType::Pointer transform) {
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	ResamplerType::Pointer resampler = ResamplerType::New();
	resampler->SetInput(src);
	resampler->SetTransform(
                            static_cast<ResamplerType::TransformType::Pointer>(transform));
	resampler->SetInterpolator(interpolator);
	resampler->SetReferenceImage(src);
	resampler->UseReferenceImageOn();
	resampler->Update();

	return resampler->GetOutput();
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

int main_test_sampling(int argc, char* argv[]) {
	itkcmds::itkImageIO<ImageType> imageIO;
	ImageType::Pointer srcImg = imageIO.ReadImageT(argv[1]);
	ImageType::Pointer dstImg = imageIO.ReadImageT(argv[2]);

	Timer timer;
	timer.start();
	computeMSE(srcImg, dstImg);
	timer.stop();

	cout << "Elapsed Time: " << timer.getElapsedTimeInMilliSec() << endl;

	timer.start();
	TransformType::Pointer transform = TransformType::New();
	transform->Rotate(0, 1, angle2rad(5));
	ImageType::Pointer movingImg1 = resampleImageITK(srcImg, transform);
	timer.stop();
	cout << "Elapsed Time: " << timer.getElapsedTimeInMilliSec() << endl;

	timer.start();
	ImageType::Pointer movingImg2 = resampleImageRaw(srcImg, srcImg, transform);
	timer.stop();
	cout << "Elapsed Time: " << timer.getElapsedTimeInMilliSec() << endl;

    return 0;
}

int main_metric_eval(int argc, char* argv[]) {
	itkcmds::itkImageIO<ImageType> imageIO;
	ImageType::Pointer srcImg = imageIO.ReadImageT(argv[1]);
	ImageType::Pointer dstImg = imageIO.ReadImageT(argv[2]);

    Timer timer;
    timer.start();
	TransformType::Pointer transform = TransformType::New();
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	MetricType::Pointer metric = MetricType::New();
	metric->SetFixedImage(dstImg);
	metric->SetFixedImageRegion(dstImg->GetLargestPossibleRegion());
	metric->SetMovingImage(srcImg);
	metric->SetInterpolator(interpolator);
	metric->SetUseAllPixels(false);
	metric->SetTransform(static_cast<MetricType::TransformPointer>(transform));
	metric->SetNumberOfThreads(8);
	metric->SetFixedImageSamplesIntensityThreshold(0.5);

	try {
		metric->Initialize();
	}
	catch (itkException &ex) {
		cout << ex << endl;
	}

    timer.stop();
    cout << "Initialization: " << timer.getElapsedTimeInMilliSec() << " ms" << endl;

	transform->Rotate(0, 1, 8 * itk::Math::pi / 180.);

	MetricType::ParametersType params = transform->GetParameters();
	MetricType::DerivativeType derivOut;

    timer.start();
	double m = metric->GetValue(params);
    timer.stop();
    cout << "GetValue(): " << timer.getElapsedTimeInMilliSec() << " ms" << endl;

    timer.start();
	metric->GetDerivative(params, derivOut);
    timer.stop();
    cout << "GetDerivative(): " << timer.getElapsedTimeInMilliSec() << " ms" << endl;


	cout << "Number of Pixels: " << metric->GetNumberOfPixelsCounted() << endl;
	cout << "Initial Params: " << params << endl;
	cout << "Metric Value: " << m << endl;
	cout << "Metric Derivative: " << derivOut << endl;

	return 0;
}

int run_optimization(ImageType::Pointer src, ImageType::Pointer dst) {
    TransformType::Pointer affine = TransformType::New();
    TransformType::ParametersType params = affine->GetParameters();

    cout << "Affine parameter: " << params << endl;
    Timer timer;
    timer.start();

    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    MetricType::Pointer metric = MetricType::New();
    metric->SetMovingImage(src);
    metric->SetFixedImage(dst);
    metric->UseAllPixelsOn();
    metric->SetInterpolator(interpolator);
    metric->SetFixedImageRegion(src->GetBufferedRegion());
    metric->SetTransform(affine);
    metric->Initialize();
    timer.stop();
    cout << "Metric Initialization Done: " << timer.getElapsedTimeInMilliSec() << " ms" << endl;

    OptimizationReporter::Pointer optiView = OptimizationReporter::New();

    MyFunction::Pointer myFun = MyFunction::New();

    OptimizerType::Pointer optimizer = OptimizerType::New();
    optimizer->SetCostFunction(static_cast<OptimizerType::CostFunctionType::Pointer>(metric));
    optimizer->SetInitialPosition(params);
    optimizer->AddObserver(StartEvent(), optiView);
    optimizer->AddObserver(IterationEvent(), optiView);
    optimizer->AddObserver(EndEvent(), optiView);

    cout << "Ready to optimize ..." << endl;
    optimizer->StartOptimization();
    cout << "Optimization finished ..." << endl;


    cout << "Current Cost: " << optimizer->GetCurrentCost() << endl;
    cout << "Current Parameter: " << optimizer->GetCurrentPosition() << endl;
    cout << "Current Iteration: " << optimizer->GetCurrentIteration() << endl;

    return 0;
}

//
//class MultipleRegionRegistration {
//private:
//    typedef vector<ImageType::Pointer> ImageVector;
//    MetaMetrics::Pointer _metric;
//    MetaMetrics::ParametersType _currentParams;
//    ImageType::Pointer _fixedImage;
//    ImageVector _movingImages;
//    itkcmds::itkImageIO<ImageType> _imageIO;
//    OptimizerType::Pointer _optimizer;
//    OptimizationReporter::Pointer _reporter;
//
//public:
//    void loadImages(int argc, char* argv[]) {
//        ImageType::Pointer fixed = _imageIO.ReadImageT(argv[1]);
//        for (int i = 2; i < argc; i++) {
//            ImageType::Pointer moving = _imageIO.ReadImageT(argv[i]);
//            _movingImages.push_back(moving);
//        }
//        _reporter = OptimizationReporter::New();
//    }
//
//    void initialize() {
//        _metric = MetaMetrics::New();
//        MetaMetrics::ParametersType initialParams;
//        initialParams.SetSize(_metric->GetNumberOfParameters());
//        int nOffset = 0;
//        for (ImageVector::size_type i = 0; i < _movingImages.size(); i++) {
//            MetricType::Pointer metric = MetricType::New();
//            InterpolatorType::Pointer interpolator = InterpolatorType::New();
//            TransformType::Pointer transform = TransformType::New();
//            metric->SetFixedImage(_fixedImage);
//            metric->SetMovingImage(_movingImages[i]);
//            metric->SetFixedImageRegion(_fixedImage->GetBufferedRegion());
//            metric->UseAllPixelsOn();
//            metric->SetInterpolator(interpolator);
//            metric->SetTransform(transform);
//            metric->Initialize();
//            _metric->AddMetric(metric);
//            int nParam = transform->GetNumberOfParameters();
//            TransformType::ParametersType subParam = transform->GetParameters();
//            for (int i = 0; i < nParam; i++) {
//                initialParams[nOffset + i] = subParam[i];
//            }
//        }
//        _optimizer->SetCostFunction(_metric);
//        _optimizer->SetInitialPosition(initialParams);
//        _optimizer->AddObserver(StartEvent(), _reporter);
//        _optimizer->AddObserver(IterationEvent(), _reporter);
//        _optimizer->AddObserver(EndEvent(), _reporter);
//    }
//
//    void startOptimization() {
//        _optimizer->StartOptimization();
//    }
//};

//void cropOrientedBoundingBox(ImageType::Pointer image, int label) {
//    LabelGeometryFilter::Pointer labelFilter = LabelGeometryFilter::New();
//    labelFilter->SetInput(image);
//    labelFilter->CalculateOrientedBoundingBoxOn();
//    labelFilter->CalculateOrientedLabelRegionsOn();
//    LabelGeometryFilter::BoundingBoxVerticesType vertices = labelFilter->GetOrientedBoundingBoxVertices(3);
//    for (LabelGeometryFilter::BoundingBoxVerticesType::size_type i = 0; i < vertices.size(); i++) {
//        cout << vertices[i] << endl;
//    }
//}

ImageType::Pointer cropBoundingBox(ImageType::Pointer intensityImage, ImageType::Pointer labelImage, int label) {
    typedef itksys::hash_map< int, int> MapType;
    MapType map;
    map[1] = 1;
    map[2] = 3;
    cout << map.size() << endl;

    ImageType::IndexType lowerIndex, upperIndex;
    ImageType::RegionType cropRegion;
    cropRegion.SetIndex(lowerIndex);
    cropRegion.SetUpperIndex(upperIndex);
    CropFilterType::Pointer cropFilter = CropFilterType::New();
    cropFilter->SetInput(intensityImage);
    cropFilter->SetExtractionRegion(cropRegion);
    cropFilter->Update();
    return cropFilter->GetOutput();
    return intensityImage;
}

//
//int run_multiple_optimization(int argc, char* argv[]) {
//    MultipleRegionRegistration registration;
//    registration.loadImages(argc, argv);
//    registration.initialize();
//    registration.startOptimization();
//    return 0;
//}
//


int main_single_frpr_optimization(int argc, char* argv[]) {
    itkcmds::itkImageIO<ImageType> imageIO;
    ImageType::Pointer src = imageIO.ReadImageT(argv[1]);
    ImageType::Pointer dst = imageIO.ReadImageT(argv[2]);

    cout << "Run FRPR Optimization ..." << endl;
    run_optimization(src, dst);
    return 0;
}

int main_function_opti(int argc, char* argv[]) {
    MyFunction::Pointer myfun = MyFunction::New();

    itk::MyFRPROptimizer::Pointer opti = itk::MyFRPROptimizer::New();
    MyFunction::ParametersType param;
    param.SetSize(2);
    param[0] = 10;
    param[1] = 5;

    opti->SetCostFunction(myfun);
    opti->SetInitialPosition(param);
    opti->StartOptimization();
    cout << "Number of Iterations: " << opti->GetCurrentIteration() << endl;
    cout << "Solution: " << opti->GetCurrentPosition() << endl;

    return 0;
}
*/
