#include "itkCastImageFilter.h"
#include "itkImage.h"
#include "itkImageRegistrationMethod.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkResampleImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkTranslationTransform.h"
#include "itkRigid2DTransform.h"
#include "vnl/vnl_matlab_write.h"

#include "piImageIO.h"
#include "piOptions.h"

const    unsigned int    Dimension = 2;
typedef  unsigned char           PixelType;

typedef itk::Image< PixelType, Dimension >  ImageType;

//  The transform that will map the fixed image into the moving image.
typedef itk::TranslationTransform< double, Dimension > TransformType;
typedef itk::Rigid2DTransform<double> RigidTransformType;

//  The metric will compare how well the two images match each other. Metric
//  types are usually parameterized by the image types as it can be seen in
//  the following type declaration.
typedef itk::MeanSquaresImageToImageMetric<ImageType, ImageType> MetricType;

//  Finally, the type of the interpolator is declared. The interpolator will
//  evaluate the intensities of the moving image at non-grid positions.
typedef itk:: LinearInterpolateImageFunction<ImageType, double> InterpolatorType;


int optimizeRegistration(ImageType::Pointer fixedImage, ImageType::Pointer movingImage, ImageType::Pointer outputImage) {
    
    //  An optimizer is required to explore the parameter space of the transform
    //  in search of optimal values of the metric.
    typedef itk::RegularStepGradientDescentOptimizer       OptimizerType;
    
    
    //  The registration method type is instantiated using the types of the
    //  fixed and moving images. This class is responsible for interconnecting
    //  all the components that we have described so far.
    typedef itk::ImageRegistrationMethod<
    ImageType,
    ImageType >    RegistrationType;
    
    // Create components
    MetricType::Pointer         metric        = MetricType::New();
    RigidTransformType::Pointer      transform     = RigidTransformType::New();
    InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
    OptimizerType::Pointer      optimizer     = OptimizerType::New();
    RegistrationType::Pointer   registration  = RegistrationType::New();
    
    // Each component is now connected to the instance of the registration method.
    registration->SetMetric(        metric        );
    registration->SetOptimizer(     optimizer     );
    registration->SetTransform(     transform     );
    registration->SetInterpolator(  interpolator  );
    
    
   // Set the registration inputs
    registration->SetFixedImage(fixedImage);
    registration->SetMovingImage(movingImage);
    registration->SetFixedImageRegion(
                                      fixedImage->GetLargestPossibleRegion() );
    
    //  Initialize the transform
    typedef RegistrationType::ParametersType ParametersType;
    ParametersType initialParameters( transform->GetNumberOfParameters() );
    
    initialParameters[0] = -5*3.141592/180.0;  // Initial offset along X
    initialParameters[1] = -40;  // Initial offset along Y
    initialParameters[2] = 0;  // Initial offset along Y
    
    registration->SetInitialTransformParameters( initialParameters );
    
    optimizer->SetMaximumStepLength( 4.00 );
    optimizer->SetMinimumStepLength( 0.01 );
    
    // Set a stopping criterion
    optimizer->SetNumberOfIterations( 200 );
    
    // Connect an observer
    //CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
    //optimizer->AddObserver( itk::IterationEvent(), observer );
    
    try
    {
        registration->Update();
    }
    catch( itk::ExceptionObject & err )
    {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
    }
    
    //  The result of the registration process is an array of parameters that
    //  defines the spatial transformation in an unique way. This final result is
    //  obtained using the \code{GetLastTransformParameters()} method.
    
    ParametersType finalParameters = registration->GetLastTransformParameters();
    
    //  In the case of the \doxygen{TranslationTransform}, there is a
    //  straightforward interpretation of the parameters.  Each element of the
    //  array corresponds to a translation along one spatial dimension.
    
    const double TranslationAlongX = finalParameters[0];
    const double TranslationAlongY = finalParameters[1];
    
    //  The optimizer can be queried for the actual number of iterations
    //  performed to reach convergence.  The \code{GetCurrentIteration()}
    //  method returns this value. A large number of iterations may be an
    //  indication that the maximum step length has been set too small, which
    //  is undesirable since it results in long computational times.
    
    const unsigned int numberOfIterations = optimizer->GetCurrentIteration();
    
    //  The value of the image metric corresponding to the last set of parameters
    //  can be obtained with the \code{GetValue()} method of the optimizer.
    
    const double bestValue = optimizer->GetValue();
    
    // Print out results
    //
    std::cout << "Result = " << std::endl;
    std::cout << " Translation X = " << TranslationAlongX  << std::endl;
    std::cout << " Translation Y = " << TranslationAlongY  << std::endl;
    std::cout << " Iterations    = " << numberOfIterations << std::endl;
    std::cout << " Metric value  = " << bestValue          << std::endl;
    
    //  It is common, as the last step of a registration task, to use the
    //  resulting transform to map the moving image into the fixed image space.
    //  This is easily done with the \doxygen{ResampleImageFilter}. Please
    //  refer to Section~\ref{sec:ResampleImageFilter} for details on the use
    //  of this filter.  First, a ResampleImageFilter type is instantiated
    //  using the image types. It is convenient to use the fixed image type as
    //  the output type since it is likely that the transformed moving image
    //  will be compared with the fixed image.
    
    typedef itk::ResampleImageFilter<
    ImageType,
    ImageType >    ResampleFilterType;
    
    //  A resampling filter is created and the moving image is connected as  its input.
    
    ResampleFilterType::Pointer resampler = ResampleFilterType::New();
    resampler->SetInput( movingImage);
    
    //  The Transform that is produced as output of the Registration method is
    //  also passed as input to the resampling filter. Note the use of the
    //  methods \code{GetOutput()} and \code{Get()}. This combination is needed
    //  here because the registration method acts as a filter whose output is a
    //  transform decorated in the form of a \doxygen{DataObject}. For details in
    //  this construction you may want to read the documentation of the
    //  \doxygen{DataObjectDecorator}.
    
    resampler->SetTransform( registration->GetOutput()->Get() );
    
    //  As described in Section \ref{sec:ResampleImageFilter}, the
    //  ResampleImageFilter requires additional parameters to be specified, in
    //  particular, the spacing, origin and size of the output image. The default
    //  pixel value is also set to a distinct gray level in order to highlight
    //  the regions that are mapped outside of the moving image.
    
    resampler->SetSize( fixedImage->GetLargestPossibleRegion().GetSize() );
    resampler->SetOutputOrigin(  fixedImage->GetOrigin() );
    resampler->SetOutputSpacing( fixedImage->GetSpacing() );
    resampler->SetOutputDirection( fixedImage->GetDirection() );
    resampler->SetDefaultPixelValue( 100 );
    
    //  The output of the filter is passed to a writer that will store the
    //  image in a file. An \doxygen{CastImageFilter} is used to convert the
    //  pixel type of the resampled image to the final type used by the
    //  writer. The cast and writer filters are instantiated below.
    
    outputImage = resampler->GetOutput();
    
    return EXIT_SUCCESS;
}


void generateMetricMap(ImageType::Pointer fixedImage, ImageType::Pointer movingImage, string outputFile) {
    MetricType::Pointer         metric        = MetricType::New();
    RigidTransformType::Pointer      transform     = RigidTransformType::New();
    InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
    
    TransformType::ParametersType params;
    params.Initialize();
    params.SetSize(3);
    
    
    RigidTransformType::ParametersType fixedParams;
    fixedParams.SetSize(2);
    fixedParams[0] = 256;
    fixedParams[1] = 256;
    transform->SetFixedParameters(fixedParams);
    
    // set up metric with given arguments
    metric->SetFixedImage(fixedImage);
    metric->SetFixedImageRegion(fixedImage->GetBufferedRegion());
    metric->SetMovingImage(movingImage);
    metric->SetTransform(transform);
    metric->SetInterpolator(interpolator);
    
    const int N = 50;
    const int D = 30;
    vnl_matrix<double> mat(2*D+1, 2*N+1);
    mat.fill(0);
    for (int j = -N; j <= N; j+=2) {
        for (double k = -D; k <= D; k+=2.5) {

            params[0] = k * 3.141592/180.0;
            params[1] = j;
            params[2] = 0;
            
            metric->SetUseSequentialSampling(true);
            metric->UseAllPixelsOn();
            metric->Initialize();
            
            double value = metric->GetValue(params);
            mat(k+D,j+N) = value;
            cout << k << "(deg) ," << j << " (px) = " << value << endl;
        }
    }

    ofstream of(outputFile);
    vnl_matlab_write(of, mat.data_array(), mat.rows(), mat.cols(), "MSE");
    of.close();
}


void applyTransform(ImageType::Pointer inputImage, ImageType::Pointer& outputImage, pi::IntVector intParams) {
    typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
    RigidTransformType::Pointer      transform     = RigidTransformType::New();
    RigidTransformType::ParametersType params;
    params.SetSize(3);
    params[0] = intParams[0]*3.141592/180.0;
    params[1] = intParams[1];
    params[2] = intParams[2];
    transform->SetParameters(params);
    
    cout << "Txf Parameters: " << params << endl;
    
    RigidTransformType::ParametersType fixedParams;
    fixedParams.SetSize(2);
    fixedParams[0] = 256;
    fixedParams[1] = 256;
    transform->SetFixedParameters(fixedParams);
    cout << "Txf Fixed Parameters: " << fixedParams << endl;

//    transform->SetFixedParameters(params);

    
    //  A resampling filter is created and the moving image is connected as  its input.
    
    ResampleFilterType::Pointer resampler = ResampleFilterType::New();
    resampler->SetInput( inputImage);
    
    //  The Transform that is produced as output of the Registration method is
    //  also passed as input to the resampling filter. Note the use of the
    //  methods \code{GetOutput()} and \code{Get()}. This combination is needed
    //  here because the registration method acts as a filter whose output is a
    //  transform decorated in the form of a \doxygen{DataObject}. For details in
    //  this construction you may want to read the documentation of the
    //  \doxygen{DataObjectDecorator}.
    
    resampler->SetTransform( transform );
    
    //  As described in Section \ref{sec:ResampleImageFilter}, the
    //  ResampleImageFilter requires additional parameters to be specified, in
    //  particular, the spacing, origin and size of the output image. The default
    //  pixel value is also set to a distinct gray level in order to highlight
    //  the regions that are mapped outside of the moving image.
    
    resampler->SetSize( inputImage->GetLargestPossibleRegion().GetSize() );
    resampler->SetOutputOrigin(  inputImage->GetOrigin() );
    resampler->SetOutputSpacing( inputImage->GetSpacing() );
    resampler->SetOutputDirection( inputImage->GetDirection() );
    resampler->SetDefaultPixelValue( 255 );
    
    //  The output of the filter is passed to a writer that will store the
    //  image in a file. An \doxygen{CastImageFilter} is used to convert the
    //  pixel type of the resampled image to the final type used by the
    //  writer. The cast and writer filters are instantiated below.
    
    resampler->Update();
    outputImage = resampler->GetOutput();
    
    cout << "Resample done..." << outputImage.GetPointer() << endl;
    
    pi::ImageIO<ImageType> imgIO;
    imgIO.WriteImage("output.png", outputImage);

    
    
}


int main(int argc, char* argv[]) {
    pi::Options opts;
    // general options
    opts.addOption("--reg", "run registration", "fixed-image moving-image output-image", SO_NONE);
    opts.addOption("--mse", "generate meansquare error map", "fixed-image moving-image output-image", SO_NONE);
    opts.addOption("--txf", "generate transformed image", "moving-image --params=angle,tx,ty", SO_NONE);
    opts.addOption("--params", "set parameter values", SO_REQ_SEP);
    
    pi::StringVector args = opts.ParseOptions(argc, argv, NULL);

    pi::ImageIO<ImageType> imgIO;
    ImageType::Pointer fixedImage = imgIO.ReadCastedImage(args[0]);
    ImageType::Pointer movingImage;
    if (args.size() > 1) {
        movingImage = imgIO.ReadCastedImage(args[1]);
    }
    ImageType::Pointer outputImage;
    
    if (opts.GetBool("--reg")) {
        optimizeRegistration(fixedImage, movingImage, outputImage);
        imgIO.WriteImage(args[2], outputImage);
    } else if (opts.GetBool("--mse")) {
        cout << "generating MSE metric map" << endl;
        generateMetricMap(fixedImage, movingImage, args[2]);
        if (!outputImage.IsNull()) {
            imgIO.WriteImage(args[2], outputImage);
        }
    } else if (opts.GetBool("--txf")) {
        pi::IntVector params = opts.GetStringAsIntVector("--params");
        applyTransform(fixedImage, outputImage, params);
        if (!outputImage.IsNull()) {
            imgIO.WriteImage(args[1], outputImage);
            cout << "writing done ..." << endl;
        }
    }
}
