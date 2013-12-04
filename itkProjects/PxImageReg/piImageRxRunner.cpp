//
//  piImageRxRunner.cpp
//  PxImageReg
//
//  Created by Joohwi Lee on 11/30/13.
//
//

#include "piImageRxRunner.h"
#include "piOptions.h"
#include "piImageIO.h"
#include "piImageDef.h"
#include "piImageProc.h"

#include <itkEuler2DTransform.h>
#include <itkEuler3DTransform.h>
#include <itkMeanSquaresImageToImageMetricv4.h>
#include <itkCorrelationImageToImageMetricv4.h>
#include <itkMattesMutualInformationImageToImageMetricv4.h>
#include <itkRegistrationParameterScalesFromIndexShift.h>
#include <itkGradientDescentOptimizerv4.h>
#include "itk/itkEntropyImageToImageMetricv4.h"

namespace pi {
    void executeRxRunner(Options& opts, StringVector& args) {
        ImageRx runner(opts, args);
        runner.mainRigidRegistration(opts, args);
        exit(0);
    }


#pragma mark MovingImage Implementations
    MovingImage::MovingImage()
    {
        _transform = NULL;
    }

    void MovingImage::setImage(RealImage::Pointer image) {
        _image = image;

        // use the half size of spacing as a sigma
        double sigma = _image->GetSpacing()[0] / 2.0;

        // compute the gradient image using ImageProc
        _gradientImage = ComputeGaussianGradient(_image, sigma);

        // set up interpolators
        _intensityInterpolator = itk::LinearInterpolateImageFunction<RealImage>::New();
        _intensityInterpolator->SetInputImage(_image);

        _gradientInterpolator = itk::VectorLinearInterpolateImageFunction<GradientImage>::New();
        _gradientInterpolator->SetInputImage(_gradientImage);
    }

    RealImage::PixelType MovingImage::samplePixel(RealImage::PointType fixedPoint, bool& isValidPixel) {
        // first compute the transformed point at this moving image
        if (_transform) {
            RealImage::PointType movingPoint = _transform->TransformPoint(fixedPoint);
            if (_intensityInterpolator->IsInsideBuffer(movingPoint)) {
                isValidPixel = true;
                return _intensityInterpolator->Evaluate(movingPoint);
            } else {
                // if the moving point is outside of the image buffer
                // the return pixel is not defined
                isValidPixel = false;
                return 0;
            }
        } else {
            return _intensityInterpolator->Evaluate(fixedPoint);
        }
    }


    GradientImage::PixelType MovingImage::sampleGradient(RealImage::PointType fixedPoint, bool& isValidPixel) {
        // first compute the transformed point at this moving image
        if (_transform) {
            RealImage::PointType movingPoint = _transform->TransformPoint(fixedPoint);
            if (_gradientInterpolator->IsInsideBuffer(movingPoint)) {
                isValidPixel = true;
                return _gradientInterpolator->Evaluate(movingPoint);
            } else {
                // if the moving point is outside of the image buffer
                // the return gradient is not defined
                isValidPixel = false;
                return GradientImage::PixelType();
            }
        } else {
            return _gradientInterpolator->Evaluate(fixedPoint);
        }
    }

    // set a transform instance
    void MovingImage::setTransform(TransformType *transform) {
        _transform = transform;
    }

    // provide parameters to the transform instance
    void MovingImage::setTransformParameters(TransformType::ParametersType params) {
        if (_transform) {
            _transform->SetParameters(params);
        }
    }


#pragma mark EntropyImageMetric Implementations
void EntropyImageMetric::setFixedImage(RealImage::Pointer image) {
    // the fixedImage will not have transform
    _fixedImage.setImage(image);
}

void EntropyImageMetric::addMovingImage(RealImage::Pointer image, TransformType* transform, TransformType::ParametersType& params) {
    // the movingImage will have its own transform and its initial parameters
    MovingImage movingImage;
    movingImage.setImage(image);
    movingImage.setTransform(transform);
    movingImage.setTransformParameters(params);

    // add to moving image vectors so that additional images can be added
    _movingImages.push_back(movingImage);
}

double EntropyImageMetric::GetValue() const {
    // get value will return \log \det (Covariance)
    return 0;
}

void EntropyImageMetric::GetValueAndDerivative(double &value, DerivativeType &deriv) const {
    // compute value and its derivatives
    return;
}


#pragma mark OptimizerProgress Listener

    /**
     * This class is called from a registration module in order to provide internal registration information to a user
     *
     */
    template <class T>
    class OptimizerProgress: public itk::Command {
    public:
        /** Standard class typedefs. */
        typedef OptimizerProgress Self;
        typedef itk::Command Superclass;
        typedef itk::SmartPointer<Self> Pointer;
        typedef itk::SmartPointer<const Self> ConstPointer;

        itkTypeMacro(OptimizerProgress, itk::Command);
        itkNewMacro(OptimizerProgress);

        virtual void Execute(Object *caller, const itk::EventObject & event) {
            T* realCaller = dynamic_cast<T*>(caller);
            if (realCaller == NULL) {
                return;
            }
            cout << realCaller->GetCurrentIteration() << ": " << realCaller->GetCurrentPosition() << endl;
            realCaller->Print(cout);
        }

        /**
         * when a caller is defined as const
         */
        virtual void Execute(const Object *caller, const itk::EventObject & event) {
        }

    protected:
        OptimizerProgress() {}
        virtual ~OptimizerProgress() {}

    private:
        // purposely not implemented
        OptimizerProgress(const Self &);
        void operator=(const Self &);
    };




#pragma mark Registration Functions

    /***
     * perform rigid registration
     *
     * command line:
     * --rx fixedImage movingImage [resampledImage] [transformParameters]
     *
     */

    void ImageRx::mainRigidRegistration(Options& opts, StringVector& args) {

        if (args.size() < 4) {
            cout << "--rx [fixed-image] [moving-image] [output-image] [output-transform]" << endl;
            return;
        }

#define USE_AFFINE_TRANSFORM

#ifdef USE_RIGID_TRANSFORM
#if DIMENSIONS == 2
        typedef itk::Euler2DTransform<> TransformType;
#else
        typedef itk::Euler3DTransform<> TransformType;
#endif
#endif

#ifdef USE_AFFINE_TRANSFORM
        typedef itk::AffineTransform<double,DIMENSIONS> TransformType;
#endif

        ImageIO<RealImage> imageIO;

        RealImage::Pointer fixedImage  = imageIO.ReadCastedImage(args[0]);
        RealImage::Pointer movingImage = imageIO.ReadCastedImage(args[1]);


        // There are three main components in the ITK registration framework.

        // 1) Optimizer computes a gradient descent direction and updates the parameters toward a minimizing direction.
        // 2) Similarity Metric provides a cost function as well as derivative functions.
        // 3) Transformation and its parameters transform the fixed image space to moving image so that a resampled image can be computed
        //

        //        typedef itk::TranslationTransform<double,3> TransformType;
        //        typedef itk::VersorRigid3DTransform<double> TransformType;
        //        typedef itk::Similarity3DTransform<double> TransformType;
        //        typedef itk::BSplineTransform<double,2,4> TransformType;
        //        typedef itk::AffineTransform<double,3> TransformType;
        TransformType::Pointer transform;
        transform = TransformType::New();

        // define a gradient descent optimizer function
        typedef itk::GradientDescentOptimizerv4 OptimizerType;
        OptimizerType::Pointer optimizer;

        // define a MeanSquares image to image metric
        typedef OptimizerType::ParametersType ParametersType;





        // set up initial parameters
        ParametersType initialParams;
        initialParams.SetSize(transform->GetNumberOfParameters());
        initialParams.Fill(0);

#ifdef USE_AFFINE_TRANSFORM
        initialParams[0] = 1;
        initialParams[3] = 1;
#endif


        // set up fixed parameters: the center of rotation
        ParametersType fixedParams;
        fixedParams.SetSize(fixedImage->GetImageDimension());

        itk::ContinuousIndex<double,RealImage::ImageDimension> centerIdx;
        RealImage::PointType centerPoint;

        for (int i = 0; i < fixedParams.GetSize(); i++) {
            centerIdx[i] = fixedImage->GetBufferedRegion().GetIndex(i) + fixedImage->GetBufferedRegion().GetSize(i) / 2.0;
        }
        fixedImage->TransformContinuousIndexToPhysicalPoint(centerIdx, centerPoint);
        for (int i = 0; i < fixedParams.size(); i++) {
            fixedParams[i] = centerPoint[i];
        }
        transform->SetFixedParameters(fixedParams);


        typedef itk::LinearInterpolateImageFunction<RealImage> ImageInterpolator;

        // set up cost functions
        // the parameters are different for each cost function

#define USE_MSQ

#ifdef USE_MSQ
        typedef itk::MeanSquaresImageToImageMetricv4<RealImage, RealImage> CostFunctionType;
#endif
#ifdef USE_CC
        typedef itk::CorrelationImageToImageMetricv4<RealImage, RealImage> CostFunctionType;
#endif


        CostFunctionType::Pointer costFunc;
        costFunc = CostFunctionType::New();
        costFunc->SetFixedImage(fixedImage);
        costFunc->SetFixedInterpolator(ImageInterpolator::New());
        costFunc->SetMovingImage(movingImage);
        costFunc->SetMovingInterpolator(ImageInterpolator::New());
        costFunc->SetMovingTransform(transform);
        costFunc->SetParameters(initialParams);
        costFunc->Initialize();


#ifdef USE_MI
        typedef itk::MattesMutualInformationImageToImageMetricv4<RealImage, RealImage> CostFunctionType;
        CostFunctionType::Pointer costFunc;
        costFunc = CostFunctionType::New();
        costFunc->SetFixedImage(fixedImage);
        costFunc->SetFixedInterpolator(ImageInterpolator::New());
        costFunc->SetMovingImage(movingImage);
        costFunc->SetMovingInterpolator(ImageInterpolator::New());
        costFunc->SetMovingTransform(transform);
        costFunc->SetParameters(initialParams);
        if (dynamic_cast<itk::MattesMutualInformationImageToImageMetricv4<RealImage, RealImage>*>(costFunc.GetPointer()) != NULL) {
            costFunc->SetNumberOfHistogramBins(32);
        }
        costFunc->Initialize();
#endif

        typedef itk::RegistrationParameterScalesFromIndexShift<CostFunctionType> ScaleEstimatorType;
        ScaleEstimatorType::Pointer estimator = ScaleEstimatorType::New();
        estimator->SetMetric(costFunc);

        // print out optimizer progress to console
        OptimizerProgress<OptimizerType>::Pointer progress = OptimizerProgress<OptimizerType>::New();

        // normalize paramter update scaling
        OptimizerType::ScalesType scales;
        scales.SetSize(transform->GetNumberOfParameters());
        scales.Fill(1);


        // set up the optimizer
        optimizer = OptimizerType::New();
        optimizer->SetScales(scales);
        optimizer->SetScalesEstimator(estimator);
        optimizer->SetMetric(costFunc);

        // number of iterations
        optimizer->SetNumberOfIterations(1000);
        RealImage::SpacingType spacing = fixedImage->GetSpacing();
        optimizer->SetMaximumStepSizeInPhysicalUnits(spacing[0]*3);

        // set the optimizer observer
        optimizer->AddObserver(itk::IterationEvent(), progress);

        try {
            ::itk::Object::GlobalWarningDisplayOn();
            optimizer->SetDebug(true);
            costFunc->SetDebug(true);
            optimizer->Print(cout);
            optimizer->StartOptimization();
        } catch (itk::ExceptionObject& ex) {
            ex.Print(cout);
        }

        cout << "Iterations: " << optimizer->GetCurrentIteration() << endl;
        cout << "Stop Reason: " << optimizer->GetStopConditionDescription() << endl;

        const ParametersType& params = optimizer->GetCurrentPosition();
        transform->SetParameters(params);

        typedef itk::ResampleImageFilter<RealImage,RealImage> ResampleFilter;
        ResampleFilter::Pointer resampler = ResampleFilter::New();
        resampler->SetInput(movingImage);
        if (transform.IsNotNull()) {
            cout << "Setting transform:" << endl;
            transform->Print(cout);
            resampler->SetTransform(dynamic_cast<ResampleFilter::TransformType*>(transform.GetPointer()));
        }
        resampler->SetReferenceImage(fixedImage);
        resampler->UseReferenceImageOn();
        resampler->Update();

        RealImage::Pointer resampled = resampler->GetOutput();
        imageIO.WriteImage(args[2], resampled);
        imageIO.WriteSingleTransform(args[3].c_str(), transform);
    }
}
