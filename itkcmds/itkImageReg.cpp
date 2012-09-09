/*
 * itkImageReg.cpp
 *
 *  Created on: Sep 5, 2012
 *      Author: joohwi
 */

#include "itkImageCommon.h"
#include "itkExceptionObject.h"
#include "itkAffineTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkImageRegistrationMethod.h"
#include "itkResampleImageFilter.h"
#include "itkContinuousIndex.h"
#include "itkMath.h"
#include "itkTimer.h"
#include "iostream"
#include "vector"
#include "itkMyMetric.h"
#include "itkMyFRPROptimizer.h"
#include "MyFunction.h"
#include "MatrixCode.h"
#include "itkCommand.h"
#include "itkSingleValuedNonLinearOptimizer.h"
#include "itkRealTimeClock.h"


using namespace std;
using namespace itk;

typedef itk::ExceptionObject itkException;
typedef itk::Image<double, 3> ImageType;
typedef itk::AffineTransform<double, 3> TransformType;
typedef itk::MyFRPROptimizer OptimizerType;
typedef itk::MyMetric<ImageType, ImageType> MetricType;
typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;
typedef itk::ImageRegistrationMethod<ImageType, ImageType> RegistrationType;
typedef itk::ResampleImageFilter<ImageType, ImageType> ResamplerType;

#define angle2rad(a) (a*itk::Math::pi/180.0)


class OptimizationReporter : public itk::Command {
public:
    typedef OptimizationReporter Self;
    typedef itk::Command Superclass;
    typedef itk::SmartPointer<Self> Pointer;

    itkTypeMacro(OptimizationReporter, itk::Command);
    itkNewMacro(OptimizationReporter);

    typedef itk::MyFRPROptimizer OptimizerType;

private:

    OptimizationReporter() {
        _updateInterval = 1;
        _numOfIterations = 0;
        _clock = itk::RealTimeClock::New();
    }

    ~OptimizationReporter() {
    }

    itk::RealTimeClock::Pointer _clock;
    itk::RealTimeClock::TimeStampType _lastTime;
    int _numOfIterations;
    int _updateInterval;

public:

    void Execute(itk::Object* caller, const itk::EventObject& event) {
        Execute((const Object*) caller, event);
    }

    void Execute(const itk::Object* object, const EventObject& event) {
        if (object == NULL) {
            cout << "Null sender is not processed..." << endl;
            return;
        }

        if (typeid(event) == typeid(itk::IterationEvent)) {
            const OptimizerType* opt = dynamic_cast<const OptimizerType*>(object);
            if (++ _numOfIterations % _updateInterval == 0) {
                itk::RealTimeClock::TimeStampType t = _clock->GetTimeInSeconds();
                cout << _numOfIterations << "\t" << opt->GetValue() << "\t" << opt->GetCurrentPosition() << "\t" << (t - _lastTime) << " secs" << endl;
                _lastTime = t;
            }
        } else if (typeid(event) == typeid(StartEvent)) {
            cout << "Optimization has started ..." << endl;
            _lastTime = _clock->GetTimeInSeconds();
        } else if (typeid(event) == typeid(EndEvent)) {
            cout << "Optimization has finisehd ..." << endl;
        }
    }

    void Update() {
        this->Execute((const Object*) NULL, IterationEvent());
    }
};

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

class MetaMetrics : public SingleValuedCostFunction {
public:
    typedef vector<MetricType::Pointer> MetricArrayType;
    typedef vector<int> IntArrayType;

private:
    MetricArrayType _metrics;
    IntArrayType _metricParamterIndex;

public:
    /** Standard class typedefs. */
    typedef MetaMetrics   Self;
    typedef SingleValuedCostFunction               Superclass;
    typedef SmartPointer< Self >       Pointer;
    typedef SmartPointer< const Self > ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(MetaMetrics, SingleValuedCostFunction);

    /**  MeasureType typedef.
     *  It defines a type used to return the cost function value. */
    typedef double MeasureType;

    /**  ParametersType typedef.
     *  It defines a position in the optimization search space. */
    typedef Superclass::ParametersType      ParametersType;
    typedef Superclass::ParametersValueType ParametersValueType;

    /** DerivativeType typedef.
     *  It defines a type used to return the cost function derivative.  */
    typedef Array< ParametersValueType > DerivativeType;

    /** Metric Vector typedef.
     *  It defines a type used to store multiple number of metrics. */
    void AddMetric(MetricType::Pointer metric) {
        _metrics.push_back(metric);
        _metricParamterIndex.push_back(metric->GetNumberOfParameters());
    }


    /** Return the number of parameters required to compute
     *  this cost function.
     *  This method MUST be overloaded by derived classes. */
    virtual unsigned int GetNumberOfParameters(void) const  {
        int nParams = 0;
        for (MetricArrayType::size_type i = 0 ; i < _metrics.size(); i++) {
            nParams += _metrics[i]->GetNumberOfParameters();
        }
        return nParams;
    }

    /** This method returns the value of the cost function corresponding
     * to the specified parameters.    */
    virtual MeasureType GetValue(const ParametersType & parameters) const {
        MeasureType value = 0;
        int nOffset = 0;
        for (MetricArrayType::size_type i = 0; i < _metrics.size(); i++) {
            int nParam = _metrics[i]->GetNumberOfParameters();
            MetricType::ParametersType param;
            param.SetSize(nParam);
            for (int j = 0; j < nParam; j++) {
                param[j] = parameters[nOffset + j];
            }
            value += _metrics[i]->GetValue(param);
            nOffset += nParam;
        }
        return value;
    }

    /** This method returns the derivative of the cost function corresponding
     * to the specified parameters.   */
    virtual void GetDerivative(const ParametersType & parameters,
                               DerivativeType & derivative) const {
        derivative.SetSize(this->GetNumberOfParameters());
        int nOffset = 0;
        for (MetricArrayType::size_type i = 0; i < _metrics.size(); i++) {
            int nParam = _metrics[i]->GetNumberOfParameters();
            MetricType::ParametersType param;
            param.SetSize(nParam);
            for (int j = 0; j < nParam; j++) {
                param[j] = parameters[nOffset + j];
            }
            MetricType::DerivativeType deriv;
            _metrics[i]->GetDerivative(param, deriv);
            for (int j = 0; j < nParam; j++) {
                derivative[nOffset + j] = deriv[j];
            }
            nOffset += nParam;
        }
    }

    /** This method returns the value and derivative of the cost function corresponding
     * to the specified parameters    */
    virtual void GetValueAndDerivative(const ParametersType & parameters,
                                       MeasureType & value,
                                       DerivativeType & derivative) const
    {
        value = this->GetValue(parameters);
        this->GetDerivative(parameters, derivative);
    }


protected:
    MetaMetrics() {}
    virtual ~MetaMetrics() {}
private:
    MetaMetrics(const Self &); //purposely not implemented
    void operator=(const Self &);           //purposely not implemented
};

class MultipleRegionRegistration {
private:
    typedef vector<ImageType::Pointer> ImageVector;
    MetaMetrics::Pointer _metric;
    MetaMetrics::ParametersType _currentParams;
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
        _metric = MetaMetrics::New();
        MetaMetrics::ParametersType initialParams;
        initialParams.SetSize(_metric->GetNumberOfParameters());
        int nOffset = 0;
        for (ImageVector::size_type i = 0; i < _movingImages.size(); i++) {
            MetricType::Pointer metric = MetricType::New();
            InterpolatorType::Pointer interpolator = InterpolatorType::New();
            TransformType::Pointer transform = TransformType::New();
            metric->SetFixedImage(_fixedImage);
            metric->SetMovingImage(_movingImages[i]);
            metric->SetFixedImageRegion(_fixedImage->GetBufferedRegion());
            metric->UseAllPixelsOn();
            metric->SetInterpolator(interpolator);
            metric->SetTransform(transform);
            metric->Initialize();
            _metric->AddMetric(metric);
            int nParam = transform->GetNumberOfParameters();
            TransformType::ParametersType subParam = transform->GetParameters();
            for (int i = 0; i < nParam; i++) {
                initialParams[nOffset + i] = subParam[i];
            }
        }
        _optimizer->SetCostFunction(_metric);
        _optimizer->SetInitialPosition(initialParams);
        _optimizer->AddObserver(StartEvent(), _reporter);
        _optimizer->AddObserver(IterationEvent(), _reporter);
        _optimizer->AddObserver(EndEvent(), _reporter);
    }

    void startOptimization() {
        _optimizer->StartOptimization();
    }
};

int run_multiple_optimization(int argc, char* argv[]) {
    MultipleRegionRegistration registration;
    registration.loadImages(argc, argv);
    registration.initialize();
    registration.startOptimization();
    return 0;
}

int main(int argc, char* argv[]) {
    run_multiple_optimization(argc, argv);
}

int run_registration(int argc, char* argv[]) {
	MetricType::Pointer metric = MetricType::New();
	TransformType::Pointer transform = TransformType::New();
	OptimizerType::Pointer optimizer = OptimizerType::New();
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	RegistrationType::Pointer registration = RegistrationType::New();

    
    itk::MyFRPROptimizer::Pointer frprOptimizer = itk::MyFRPROptimizer::New();
	metric->SetUseAllPixels(true);

	registration->SetGlobalWarningDisplay(true);
	registration->SetDebug(true);
	registration->SetMetric(metric);
	registration->SetOptimizer(frprOptimizer);
	registration->SetTransform(
			static_cast<RegistrationType::TransformPointer>(transform));
	registration->SetInterpolator(interpolator);

	itkcmds::itkImageIO<ImageType> imageIO;
	ImageType::Pointer srcImg = imageIO.ReadImageT(argv[1]);
	ImageType::Pointer dstImg = imageIO.ReadImageT(argv[2]);

	registration->SetFixedImage(dstImg);
	registration->SetMovingImage(srcImg);
	registration->SetFixedImageRegion(dstImg->GetLargestPossibleRegion());

	transform->Rotate(0, 1, 5 * itk::Math::pi / 180.);
	TransformType::ParametersType txParam = transform->GetParameters();

	RegistrationType::ParametersType initialParams(txParam);
	registration->SetInitialTransformParameters(txParam);

    /*
	optimizer->SetMaximumStepLength(4.0);
	optimizer->SetMinimumStepLength(0.01);
     */
    
	try {
		registration->Update();
		std::cout << "Registration updated ..." << std::endl;
	}
	catch (itkException &ex) {
		std::cerr << ex;
		std::cerr << std::endl;
		return EXIT_FAILURE;
	}

	return 0;
}

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