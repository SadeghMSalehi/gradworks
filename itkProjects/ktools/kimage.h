//
//  kimage.h
//  ktools
//
//  Created by Joohwi Lee on 1/22/14.
//
//

#ifndef __ktools__kimage__
#define __ktools__kimage__

#include <iostream>
#include <string>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkImageData.h>
#include <vtkDataArray.h>

#include "piImageIO.h"

typedef itk::VectorImage<float> VectorImageType;
typedef VectorImageType::PixelType VectorType;

void SetArrayTuple(vtkDataArray* a, int i, VectorType v);
void SetArrayTuple(vtkDataArray* a, int i, double v);

template <class X>
void ConvertImageT(std::string& imageFile, vtkImageData* imgData, const char* attrName, int numberOfComponents) {
    pi::ImageIO<X> itkIO;
    typename X::Pointer srcImg = itkIO.ReadCastedImage(imageFile.c_str());
    typename X::SizeType srcSize = srcImg->GetRequestedRegion().GetSize();
    typename X::PointType srcOrigin = srcImg->GetOrigin();
    typename X::SpacingType srcSpacing = srcImg->GetSpacing();
    
    imgData->SetOrigin(srcOrigin[0], srcOrigin[1], srcOrigin[2]);
    imgData->SetSpacing(srcSpacing[0], srcSpacing[1], srcSpacing[2]);
    imgData->SetDimensions(srcSize[0], srcSize[1], srcSize[2]);
    
    const int nPoints = srcImg->GetRequestedRegion().GetNumberOfPixels();
    
    vtkDoubleArray* attr = vtkDoubleArray::New();
    attr->SetNumberOfComponents(numberOfComponents);
    attr->SetName(attrName);
    attr->SetNumberOfTuples(nPoints);
    
    vtkPointData* pdata = imgData->GetPointData();
    
    switch (numberOfComponents) {
        case 1:
            pdata->SetScalars(attr);
            break;
        case 3:
            pdata->SetVectors(attr);
            break;
        case 9:
            pdata->SetTensors(attr);
            break;
        default:
            pdata->AddArray(attr);
            break;
    }
    
    int cnt = 0;
#pragma omp parallel for
    for (unsigned int z = 0; z < srcSize[2]; z++) {
        for (unsigned int y = 0; y < srcSize[1]; y++) {
            for (unsigned int x = 0; x < srcSize[0]; x++) {
                typename X::IndexType idx;
                idx[0] = x;
                idx[1] = y;
                idx[2] = z;
                typename X::PixelType v = srcImg->GetPixel(idx);
                SetArrayTuple(attr, cnt, v);
                cnt ++;
            }
        }
    }
}

#endif /* defined(__ktools__kimage__) */
