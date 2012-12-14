//
//  myImageTransform.h
//  ParticlesGUI
//
//  Created by Joohwi Lee on 12/3/12.
//
//

#ifndef __ParticlesGUI__myImageTransform__
#define __ParticlesGUI__myImageTransform__

#include <iostream>

#include "myImageContainer.h"

#include "itkTransform.h"
#include "itkKernelTransform.h"
#include "itkBSplineDeformableTransform.h"

namespace my {
    typedef itk::Transform<double, 2, 2> Transform2DType;
    typedef itk::KernelTransform<double, 2> KernelTransformType;
    typedef KernelTransformType::Pointer KernelTransformPointer;
    typedef itk::BSplineDeformableTransform<double,2,3> BSplineTransform;
    
    class ImageTransform {
    public:
        ImageTransform();
        ~ImageTransform();

        KernelTransformPointer CreateKernelTransform(int type, int n, double* src, double* dst);
        BSplineTransform::Pointer CreateBSplineTransform(int type, int n, double *src, double* dst, SliceType::Pointer spaceImage);

        //SliceType::Pointer ResampleSlice(SliceType::Pointer image, KernelTransformPointer transform);
    private:
        
    };
}

#endif /* defined(__ParticlesGUI__myImageTransform__) */
