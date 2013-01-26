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

#include "vector"

#define for4(i) for (int i = 0; i < 4; i++)
#define for3(i) for (int i = 0; i < 3; i++)

namespace my {
    const static int Dim = 3;
    
    // type definitions
    typedef itk::Image<double,Dim> DoubleImage;
    typedef itk::Image<unsigned short,Dim> LabelImage;
    typedef itk::Vector<double,Dim> VectorType;
    typedef itk::Image<VectorType,Dim> VectorImage;
    typedef itk::Offset<Dim> OffsetType;
    typedef itk::Image<OffsetType,Dim> OffsetImage;
    typedef std::vector<LabelImage::Pointer> LabelVectors;
    typedef std::vector<OffsetImage::Pointer> OffsetImageVectors;

    typedef itk::LinearInterpolateImageFunction<DoubleImage> LinearImageInterpolatorType;
    typedef itk::NearestNeighborInterpolateImageFunction<DoubleImage> NNImageInterpolatorType;
    typedef itk::VectorLinearInterpolateImageFunction<VectorImage> LinearVectorImageInterpolatorType;
    typedef itk::ImageRegionIteratorWithIndex<LabelImage> LabelImageIteratorType;


    // definition for transforms
    typedef itk::Transform<double,Dim,Dim> TransformType;

    // definition for displacement field
    typedef itk::PointSet<int,Dim> IntPointSetType;
    typedef itk::PointSet<VectorType,Dim> DisplacementFieldPointSetType;
    typedef itk::DisplacementFieldTransform<double,Dim> FieldTransformType;
    typedef itk::Image<FieldTransformType::OutputVectorType,Dim> DisplacementFieldType;

}

#endif