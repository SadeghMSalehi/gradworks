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
#include "itkMyFRPROptimizer.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkSphereBoundedGradientDescentOptimizer.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkVectorMagnitudeImageFilter.h"
#include "itkARGBColormapFunction.h"
#include "itkLBFGSBOptimizer.h"
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
typedef itk::LBFGSBOptimizer OptimizerType;
//typedef itk::SphereBoundedGradientDescentOptimizer OptimizerType;
//typedef itk::MyFRPROptimizer OptimizerType;

typedef arma::mat MatrixType;

#define USE_LBFGS_OPTIMIZER


#endif
