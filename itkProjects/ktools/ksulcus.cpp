//
//  ksulcus.cpp
//  ktools
//
//  Created by Joohwi Lee on 8/19/13.
//
//

#include "ksulcus.h"
#include "piOptions.h"
#include "piImageIO.h"
#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "itkLineSpatialObject.h"
#include "itkSpatialObjectToImageFilter.h"
#include <itkImageFileReader.h>

typedef itk::Image<unsigned char, 3> ImageType;
typedef itk::ImageFileReader<ImageType> ImageReader;
typedef itk::ImageFileWriter<ImageType> ImageWriter;

using namespace std;

ImageType::Pointer readImage(const char* filename) {
    ImageReader::Pointer reader = ImageReader::New();
    reader->SetFileName(filename);
    reader->Update();
    return reader->GetOutput();
}

void writeImage(const char* filename, ImageType::Pointer image) {
    ImageWriter::Pointer writer = ImageWriter::New();
    writer->SetFileName(filename);
    writer->SetInput(image);
    writer->Write();
    return;
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        cout << "Usage: ksulcus [vtk-model] [ref-image] [output-image]" << endl;
        return 0;
    }

    ImageType::Pointer refImage = readImage(argv[2]);
    refImage->FillBuffer(0);

    vtkPolyDataReader* reader = vtkPolyDataReader::New();
    reader->SetFileName(argv[1]);
    reader->Update();
    
    vtkPolyData* input = reader->GetOutput();

    const int nPoints = input->GetNumberOfPoints();
    for (int i = 0; i < nPoints; i++) {
        ImageType::PointType point;
        ImageType::IndexType index;

        double x[3];
        input->GetPoint(i, x);

        for (int j = 0; j < 3; j++) {
            point[j] = x[j];
        }
        refImage->TransformPhysicalPointToIndex(point, index);
        refImage->SetPixel(index, 255);
    }

    writeImage(argv[3], refImage);
}