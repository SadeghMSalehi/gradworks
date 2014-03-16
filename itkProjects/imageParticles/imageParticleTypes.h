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
#include "itkLBFGSOptimizer.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkMyPowellOptimizer.h"
#include <itkCannyEdgeDetectionImageFilter.h>
#include <itkZeroCrossingImageFilter.h>

#include "itkContinuousIndex.h"
#include "armadillo"

#define POINT_DIMENSIONS 2

const int VDimension = POINT_DIMENSIONS;
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

typedef itk::SingleValuedNonLinearOptimizer OptimizerType;
typedef OptimizerType::ParametersType OptimizerParameters;
typedef itk::RegularStepGradientDescentOptimizer GDOptimizerType;
typedef itk::LBFGSOptimizer LBFGSOptimizerType;
typedef itk::MyFRPROptimizer FRPROptimizerType;
//typedef itk::SphereBoundedGradientDescentOptimizer OptimizerType;
//typedef itk::MyPowellOptimizer OptimizerType;


typedef std::vector<float> PointVectorType;
typedef std::vector<PointVectorType> ListOfPointVectorType;
typedef std::vector<OptimizerParameters> ListOfParametersType;

typedef arma::mat MatrixType;



#define __CLAMP(x,s)((x<-3*si?0:(x>3*si?0:x)))

//#define USE_LBFGS_OPTIMIZER
//#define USE_GD_OPTIMIZER
#define USE_FRPR_OPTIMIZER
//#define USE_POWELL_OPTIMIZER



// utility functions defined in imageParticleTools.cpp

void ConvertParametersToListOfPointVectors(OptimizerParameters& inputParams, int nSubj, int nVars, ListOfPointVectorType& outputPoints);
void ConvertListOfPointVectorsToParameters(ListOfPointVectorType& inputPoints, OptimizerParameters& outputParams);
void SaveListOfPointVectors(ListOfPointVectorType& inputPoints, const char* fileName);
void LoadListOfPointVectors(const char* fileName, ListOfPointVectorType& outputPoints);



#endif
