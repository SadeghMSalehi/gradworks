//
//  myImageDef.h
//  ParticlesGUI
//
//  Created by Joohwi Lee on 1/24/13.
//
//

#ifndef ParticlesGUI_myImageDef_h
#define ParticlesGUI_myImageDef_h

#include "itkImage.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkVectorNearestNeighborInterpolateImageFunction.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkDisplacementFieldTransform.h"
#include "itkAffineTransform.h"
#include "itkCompositeTransform.h"
#include "itkTransform.h"
#include "itkGradientImageFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkPointSet.h"
#include "vnl/vnl_vector.h"
#include "piMacros.h"

#include "vector"


const static int __Dim = DIMENSIONS;

namespace pi {
    // type definitions
    typedef unsigned char LabelPixel;
    typedef float ImageReal;
    typedef double PointReal;

    typedef itk::Image<ImageReal,3> RealImage3;
    typedef itk::Image<ImageReal,2> RealImage2;
    typedef itk::Image<LabelPixel,2> LabelImage2;

    typedef itk::Image<ImageReal,__Dim> RealImage;
    typedef itk::Image<LabelPixel,__Dim> LabelImage;
    typedef itk::Vector<ImageReal,__Dim> VectorType;
    typedef itk::Image<VectorType,__Dim> VectorImage;
    typedef itk::Offset<__Dim> OffsetType;
    typedef itk::Image<OffsetType,__Dim> OffsetImage;

    typedef std::vector<OffsetImage::Pointer> OffsetImageVector;
    typedef std::vector<RealImage::Pointer> RealImageVector;
    typedef std::vector<LabelImage::Pointer> LabelImageVector;

    typedef std::vector<VectorImage::Pointer> VectorImageVector;

    typedef itk::LinearInterpolateImageFunction<RealImage> LinearImageInterpolatorType;
    typedef itk::NearestNeighborInterpolateImageFunction<RealImage> NNImageInterpolatorType;
    typedef itk::NearestNeighborInterpolateImageFunction<LabelImage> NNLabelInterpolatorType;
    typedef itk::VectorLinearInterpolateImageFunction<VectorImage> LinearVectorImageInterpolatorType;
    typedef itk::VectorNearestNeighborInterpolateImageFunction<VectorImage> NNVectorImageInterpolatorType;

    typedef NNLabelInterpolatorType::IndexType IntIndex;
    typedef LinearVectorImageInterpolatorType::ContinuousIndexType RealIndex;
    typedef RealImage::PointType ImagePoint;
    
    typedef itk::ImageRegionIteratorWithIndex<LabelImage> LabelImageIteratorType;
    typedef itk::ImageRegionIteratorWithIndex<RealImage> RealImageIteratorType;

    // gradient computation
    typedef itk::GradientImageFilter<RealImage> GradientFilterType;
    typedef GradientFilterType::OutputImageType GradientImage;
    typedef GradientFilterType::OutputPixelType GradientPixel;
    typedef itk::ImageRegionConstIteratorWithIndex<GradientImage> GradientIteratorType;
    typedef itk::ConstNeighborhoodIterator<GradientImage> VectorImageNeighborhoodIteratorType;
    typedef itk::GradientRecursiveGaussianImageFilter<RealImage, GradientImage> GaussianGradientFilterType;
    typedef itk::VectorLinearInterpolateImageFunction<GradientImage> GradientInterpolatorType;
    typedef itk::ConstNeighborhoodIterator<RealImage> RealImageNeighborhoodIteratorType;
    typedef itk::ConstNeighborhoodIterator<GradientImage> GradientImageNeighborhoodIteratorType;

    typedef itk::GradientRecursiveGaussianImageFilter<RealImage2> Gradient2FilterType;
    typedef Gradient2FilterType::OutputImageType GradientImage2;
    typedef Gradient2FilterType::OutputPixelType GradientPixel2;
    typedef itk::VectorLinearInterpolateImageFunction<GradientImage2> Gradient2InterpolatorType;
    
    // definition for transforms
    typedef itk::Transform<PointReal,__Dim,__Dim> TransformType;
    typedef itk::CompositeTransform<PointReal,__Dim> CompositeTransformType;
    typedef itk::AffineTransform<PointReal,__Dim> AffineTransformType;

    // definition for displacement field
    typedef itk::PointSet<int,__Dim> IntPointSetType;
    typedef itk::PointSet<VectorType,__Dim> DisplacementFieldPointSetType;
    typedef itk::DisplacementFieldTransform<PointReal,__Dim> FieldTransformType;
    typedef itk::Image<FieldTransformType::OutputVectorType,__Dim> DisplacementFieldType;

    // auxiliary data structures
    typedef std::vector<std::string> StringVector;


    typedef vnl_vector<double> VNLDoubleVector;
    typedef vnl_matrix<double> VNLDoubleMatrix;
    

    typedef std::vector<RealImage2::Pointer> RealImage2Vector;
    typedef std::vector<RealImage3::Pointer> RealImage3Vector;
    typedef std::vector<GradientImage2::Pointer> GradientImage2Vector;
    
    class RealImageTools {
    public:
        RealImage3::Pointer normalizeIntensity(RealImage3::Pointer image, double percentile = 0.01);
        RealImage2Vector sliceVolume(RealImage3::Pointer image, int dir);
        RealImage2Vector computeGaussianSmoothing(RealImage2Vector images, double sigma = 1);
        GradientImage2Vector computeGaussianGradient(RealImage2Vector images, double sigma = 1);
    };

    class LabelImageTools {
    public:
        std::vector<float> computeBoundingBox(LabelImage::Pointer label, int labelId);
    };

    extern RealImageTools __realImageTools;
    extern LabelImageTools __labelImageTools;
}

#endif
