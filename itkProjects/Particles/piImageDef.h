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
#include "itkPointSet.h"
#include "vnl/vnl_vector.h"

#include "vector"

#ifdef DIMENSION3
const static int __Dim = 3;
#else
const static int __Dim = 2;
#endif

#define for4(i) for (int i = 0; i < 4; i++)
#define fordim(i) for (int i = 0; i < __Dim; i++)
#define formin(x,y,z) for(int _=0;_<__Dim;_++) {z[_]=(x[_]<y[_])?x[_]:y[_];}
#define formax(x,y,z) for(int _=0;_<__Dim;_++) {z[_]=(x[_]>y[_])?x[_]:y[_];}
#define forset(x, y) for (int _=0;_<__Dim;_++) {y[_]=x[_];}
#define forcopy(x, y) for (int _=0;_<__Dim;_++) {y[_]=x[_];}
#define forfill(x, v) for (int _=0;_<__Dim;_++) {x[_]=v;}
#define forroundset(x,y) for(int _=0; _<__Dim; _++) { y[_]=x[_]+0.5; }
#define arrayset2(a,x,y) a[0]=x;a[1]=y
#define arrayset3(a,x,y,z) a[0]=x;a[1]=y;a[2]=z
#define x2string(x) x[0]<<","<<x[1]<<","<<x[2]

#ifdef DIMENSION3
#define dimdot(x,y) (x[0]*y[0]+x[1]*y[1]+x[2]*y[2])
#define dimequal(x,x0,x1,x2) (x[0]==(x0)&&x[1]==(x1)&&x[2]==(x2))
#define dimnorm2(x) (x[0]*x[0]+x[1]*x[1]+x[2]*x[2])
#else
#define dimdot(x,y) (x[0]*y[0]+x[1]*y[1])
#define dimequal(x,x0,x1) (x[0]==(x0)&&x[1]==(x1))
#define dimnorm2(x) (x[0]*x[0]+x[1]*x[1])
#endif

namespace pi {
    // type definitions
    typedef unsigned short LabelPixel;
    typedef itk::Image<double,__Dim> DoubleImage;
    typedef itk::Image<LabelPixel,__Dim> LabelImage;
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
    typedef itk::VectorNearestNeighborInterpolateImageFunction<VectorImage> NNVectorImageInterpolatorType;

    typedef NNLabelInterpolatorType::IndexType IntIndex;
    typedef LinearVectorImageInterpolatorType::ContinuousIndexType RealIndex;
    
    typedef itk::ImageRegionIteratorWithIndex<LabelImage> LabelImageIteratorType;

    // definition for transforms
    typedef itk::Transform<double,__Dim,__Dim> TransformType;
    typedef itk::CompositeTransform<double,__Dim> CompositeTransformType;
    typedef itk::AffineTransform<double,__Dim> AffineTransformType;

    // definition for displacement field
    typedef itk::PointSet<int,__Dim> IntPointSetType;
    typedef itk::PointSet<VectorType,__Dim> DisplacementFieldPointSetType;
    typedef itk::DisplacementFieldTransform<double,__Dim> FieldTransformType;
    typedef itk::Image<FieldTransformType::OutputVectorType,__Dim> DisplacementFieldType;

    // auxiliary data structures
    typedef std::vector<std::string> StringVector;

    // VNL related types
    typedef vnl_vector<double> VNLVector;
    typedef vnl_matrix<double> VNLMatrix;
}

#endif
