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
#include "itkImageRegionIteratorWithIndex.h"
#include "itkDisplacementFieldTransform.h"
#include "itkTransform.h"
#include "itkPointSet.h"
#include "vnl/vnl_vector.h"

#include "vector"

const static int __Dim = 3;

#define for4(i) for (int i = 0; i < 4; i++)
#define fordim(i) for (int i = 0; i < __Dim; i++)

namespace pi {
    
    // type definitions
    typedef itk::Image<double,__Dim> DoubleImage;
    typedef itk::Image<unsigned short,__Dim> LabelImage;
    typedef itk::Vector<double,__Dim> VectorType;
    typedef itk::Image<VectorType,__Dim> VectorImage;
    typedef itk::Offset<__Dim> OffsetType;
    typedef itk::Image<OffsetType,__Dim> OffsetImage;
    typedef std::vector<LabelImage::Pointer> LabelVector;
    typedef std::vector<OffsetImage::Pointer> OffsetImageVector;
    typedef std::vector<DoubleImage::Pointer> DoubleImageVector;
    typedef std::vector<VectorImage::Pointer> VectorImageVector;

    typedef itk::LinearInterpolateImageFunction<DoubleImage> LinearImageInterpolatorType;
    typedef itk::NearestNeighborInterpolateImageFunction<DoubleImage> NNImageInterpolatorType;
    typedef itk::NearestNeighborInterpolateImageFunction<LabelImage> NNLabelInterpolatorType;
    typedef itk::VectorLinearInterpolateImageFunction<VectorImage> LinearVectorImageInterpolatorType;
    typedef itk::ImageRegionIteratorWithIndex<LabelImage> LabelImageIteratorType;


    // definition for transforms
    typedef itk::Transform<double,__Dim,__Dim> TransformType;

    // definition for displacement field
    typedef itk::PointSet<int,__Dim> IntPointSetType;
    typedef itk::PointSet<VectorType,__Dim> DisplacementFieldPointSetType;
    typedef itk::DisplacementFieldTransform<double,__Dim> FieldTransformType;
    typedef itk::Image<FieldTransformType::OutputVectorType,__Dim> DisplacementFieldType;

    typedef std::vector<std::string> StringVector;

    // VNL related types
    typedef vnl_vector<double> VNLVector;
}

#endif
