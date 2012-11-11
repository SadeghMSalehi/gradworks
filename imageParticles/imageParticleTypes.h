//
//  imageParticleTypes.h
//  imageParticles
//
//  Created by Joohwi Lee on 10/24/12.
//
//

#ifndef imageParticles_imageParticleTypes_h
#define imageParticles_imageParticleTypes_h
#include "itkImageIO.h"
#include "itkScalarToARGBColormapImageFilter.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkMyFRPROptimizer.h"
#include "itkMyRegularStepGradientDescentOptimizer.h"
#include "itkSphereBoundedGradientDescentOptimizer.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkVectorMagnitudeImageFilter.h"
#include "itkSignedDanielssonDistanceMapImageFilter.h"
#include "itkARGBColormapFunction.h"
#include "itkLBFGSBOptimizer.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkMyPowellOptimizer.h"
#include <itkCannyEdgeDetectionImageFilter.h>
#include <itkZeroCrossingImageFilter.h>

#include "itkContinuousIndex.h"
#include "armadillo"

const int VDimension = 2;
typedef itk::CovariantVector<double,VDimension> GradientType;
typedef itk::RGBAPixel<unsigned char> RGBAPixel;
typedef itk::RGBPixel<unsigned char> RGBPixel;
typedef itk::Image<double, VDimension> ImageType;
typedef itk::Image<RGBAPixel,VDimension> BitmapType;
typedef itk::Image<GradientType,VDimension> GradientImageType;
typedef itk::GradientRecursiveGaussianImageFilter<ImageType,GradientImageType> GradientImageFilter;
typedef itk::VectorMagnitudeImageFilter<GradientImageType,ImageType> VectorMagnitudeImageFilter;
typedef itk::ScalarToARGBColormapImageFilter<ImageType, BitmapType> ScalarToRGBFilter;
typedef itk::SignedDanielssonDistanceMapImageFilter<ImageType, ImageType> DistanceMapFilter;
typedef DistanceMapFilter::VectorImageType DistanceVectorImageType;
typedef itk::LinearInterpolateImageFunction<ImageType,float> InterpolatorType;
typedef InterpolatorType::ContinuousIndexType ContinuousIndexType;
typedef itk::VectorLinearInterpolateImageFunction<DistanceVectorImageType> VectorInterpolatorType;
typedef itk::ZeroCrossingImageFilter<ImageType,ImageType> EdgeDetectionFilterType;

typedef itk::NearestNeighborInterpolateImageFunction<ImageType,float> NearestNeighborInterpolatorType;

//typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
//typedef itk::LBFGSBOptimizer OptimizerType;
//typedef itk::SphereBoundedGradientDescentOptimizer OptimizerType;
typedef itk::MyFRPROptimizer OptimizerType;
//typedef itk::MyPowellOptimizer OptimizerType;


typedef std::vector<float> PointVectorType;
typedef std::vector<PointVectorType> ListOfPointVectorType;
typedef arma::mat MatrixType;

#define __CLAMP(x,s)((x<-3*si?0:(x>3*si?0:x)))

//#define USE_LBFGS_OPTIMIZER
//#define USE_GD_OPTIMIZER
#define USE_FRPR_OPTIMIZER
//#define USE_POWELL_OPTIMIZER

#endif
