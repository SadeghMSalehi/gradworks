//
//  guiregTest.cpp
//  itktools
//
//  Created by Joohwi Lee on 9/25/12.
//
//

#include "guiregTest.h"
#include "itkImageIO.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"
#include "itkScaleVersor3DTransform.h"

void test(int argc, char* argv[]) {
    typedef itk::ScaleVersor3DTransform<double> TransformType;

    TransformType::Pointer transform = TransformType::New();
    TransformType::ParametersType param;
    param.SetSize(9);
    for (int i = 0; i < 9; i++) {
        param[i] = i;
    }
    transform->SetParameters(param);



    
}