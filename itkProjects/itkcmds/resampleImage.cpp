//
//  resampleImage.cpp
//  itkcmds
//
//  Created by Joohwi Lee on 9/19/12.
//
//

#include "itkImageIO.h"
#include "itkTransformFileReader.h"
#include "itkResampler.h"

typedef itk::Image<unsigned short,3> ImageType;
typedef itkcmds::itkImageIO<ImageType> ImageIO;

int main(int argc, char* argv[]) {
    char* outputFile = NULL;
    char* transformFile = NULL;
    ImageIO io;
    ImageType::Pointer src = io.ReadImageT(argv[1]);
    transformFile = argv[2];
    outputFile = argv[3];

    ImageType::Pointer label;
    if (argc > 4 && io.FileExists(argv[4])) {
        label = io.ReadImageT(argv[4]);
    }

    itkResampler<ImageType, ImageType, ImageIO::TransformType> resampler;
    resampler.SetInput(src);
    if (label.IsNotNull()) {
        resampler.SetInputMask(label);
    }

    ImageIO::TransformType::Pointer transform = io.ReadTransform(transformFile);
    resampler.SetTransform(transform);

    resampler.Resample();
    ImageType::Pointer resampledImage = resampler.GetOutput();

    io.WriteImageT(outputFile, resampledImage);
}