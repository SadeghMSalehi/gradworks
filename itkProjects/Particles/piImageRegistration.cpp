//
//  piImageRegistration.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/18/13.
//
//

#include "piImageRegistration.h"


#include "itkMeanSquaresImageToImageMetric.h"
#include "itkTimeProbesCollectorBase.h"
#include "itkSpatialObjectToImageFilter.h"
#include "itkEllipseSpatialObject.h"
#include "itkBSplineTransform.h"
#include "itkLBFGSOptimizer.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkSquaredDifferenceImageFilter.h"

#include <itkImageRegistrationMethod.h>
#include <itkImageRegistrationMethodv4.h>
#include <itkGradientDescentLineSearchOptimizerv4.h>


namespace pi {
    const    unsigned int    ImageDimension = DIMENSIONS;
    typedef  RealImage::PixelType           PixelType;

    typedef itk::Image< PixelType, ImageDimension >  ImageType;

    RealImage::Pointer bsplineRegistration(RealImage::Pointer srcImg, RealImage::Pointer dstImg) {

        const unsigned int SpaceDimension = ImageDimension;
        const unsigned int SplineOrder = 3;
        typedef double CoordinateRepType;

        typedef itk::BSplineTransform<CoordinateRepType, SpaceDimension, SplineOrder> TransformType;
        typedef itk::LBFGSOptimizer OptimizerType;
        typedef itk::MeanSquaresImageToImageMetric<ImageType, ImageType> MetricType;
        typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;
        typedef itk::ImageRegistrationMethod<ImageType, ImageType> RegistrationType;

        MetricType::Pointer         metric        = MetricType::New();
        OptimizerType::Pointer      optimizer     = OptimizerType::New();
        InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
        RegistrationType::Pointer   registration  = RegistrationType::New();



        // The old registration framework has problems with multi-threading
        // For now, we set the number of threads to 1
//        registration->SetNumberOfThreads(1);
        registration->SetMetric(        metric        );
        registration->SetOptimizer(     optimizer     );
        registration->SetInterpolator(  interpolator  );

        TransformType::Pointer  transform = TransformType::New();
        registration->SetTransform( transform );

        // Setup the registration
        registration->SetFixedImage(  dstImg   );
        registration->SetMovingImage(   srcImg);

        ImageType::RegionType fixedRegion = srcImg->GetBufferedRegion();
        registration->SetFixedImageRegion( fixedRegion );

        //  Here we define the parameters of the BSplineDeformableTransform grid.  We
        //  arbitrarily decide to use a grid with $5 \times 5$ nodes within the image.
        //  The reader should note that the BSpline computation requires a
        //  finite support region ( 1 grid node at the lower borders and 2
        //  grid nodes at upper borders). Therefore in this example, we set
        //  the grid size to be $8 \times 8$ and place the grid origin such that
        //  grid node (1,1) coincides with the first pixel in the fixed image.

        TransformType::PhysicalDimensionsType   fixedPhysicalDimensions;
        TransformType::MeshSizeType             meshSize;
        for (unsigned int i=0; i < ImageDimension; i++) {
            fixedPhysicalDimensions[i] = dstImg->GetSpacing()[i] *
            static_cast<double>(dstImg->GetLargestPossibleRegion().GetSize()[i] - 1 );
            meshSize[i] = dstImg->GetLargestPossibleRegion().GetSize()[i] / 8 - SplineOrder;
        }
//        unsigned int numberOfGridNodesInOneDimension = 15;
//        meshSize.Fill( numberOfGridNodesInOneDimension - SplineOrder );
        transform->SetTransformDomainOrigin( dstImg->GetOrigin() );
        transform->SetTransformDomainPhysicalDimensions( fixedPhysicalDimensions );
        transform->SetTransformDomainMeshSize( meshSize );
        transform->SetTransformDomainDirection( dstImg->GetDirection() );

        typedef TransformType::ParametersType     ParametersType;

        const unsigned int numberOfParameters = transform->GetNumberOfParameters();

        ParametersType parameters( numberOfParameters );
        parameters.Fill( 0.0 );

        transform->SetParameters( parameters );

        //  We now pass the parameters of the current transform as the initial
        //  parameters to be used when the registration process starts.

        registration->SetInitialTransformParameters( transform->GetParameters() );

        std::cout << "Intial Parameters = " << std::endl;
        std::cout << transform->GetParameters() << std::endl;

        //  Next we set the parameters of the LBFGS Optimizer.
        optimizer->SetGradientConvergenceTolerance(0.1);
        optimizer->SetLineSearchAccuracy(0.09);
        optimizer->SetDefaultStepLength(.1);
        optimizer->TraceOn();
        optimizer->SetMaximumNumberOfFunctionEvaluations(1000);

        std::cout << std::endl << "Starting Registration" << std::endl;

        try {
            registration->Update();
            std::cout << "Optimizer stop condition = "
            << registration->GetOptimizer()->GetStopConditionDescription()
            << std::endl;
        } catch (itk::ExceptionObject & err) {
            std::cerr << "ExceptionObject caught !" << std::endl;
            std::cerr << err << std::endl;
            return RealImage::Pointer();
        }

        OptimizerType::ParametersType finalParameters =
        registration->GetLastTransformParameters();
        
        std::cout << "Last Transform Parameters" << std::endl;
        std::cout << finalParameters << std::endl;
        
        transform->SetParameters( finalParameters );
        
        typedef itk::ResampleImageFilter<ImageType, ImageType>    ResampleFilterType;
        
        ResampleFilterType::Pointer resample = ResampleFilterType::New();
        
        resample->SetTransform( transform );
        resample->SetInput( srcImg );
        
        resample->SetSize(    dstImg->GetLargestPossibleRegion().GetSize() );
        resample->SetOutputOrigin(  dstImg->GetOrigin() );
        resample->SetOutputSpacing( dstImg->GetSpacing() );
        resample->SetOutputDirection( dstImg->GetDirection() );
        resample->SetDefaultPixelValue( 100 );
        resample->Update();
        return resample->GetOutput();
    }
}