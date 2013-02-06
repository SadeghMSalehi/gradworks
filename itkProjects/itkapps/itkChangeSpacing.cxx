#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

int main(int argc, char* argv[]) {
  if (argc < 5) {
    printf("usage: %s in-nrrd out-nrrd x-spacing y-spacing z-spacing [1:reset origin]\n", argv[0]);
    return 1;
  }

  bool boolResetOrigin = false;
  if (argc > 6) {
    if (1 == atoi(argv[6])) {
      boolResetOrigin = true;
    }
  }

  typedef itk::Image<double,3> ImageType;
  typedef itk::ImageFileReader<ImageType> ImageReaderType;
  typedef itk::ImageFileWriter<ImageType> ImageWriterType;

  ImageReaderType::Pointer reader = ImageReaderType::New();
  reader->SetFileName(argv[1]);
  reader->Update();

  ImageType::Pointer img = reader->GetOutput();
  double spacing[3];
  for (int d = 0; d < 3; d++) {
    spacing[d] = atof(argv[d + 3]);
  }
  img->SetSpacing(spacing);
  if (boolResetOrigin) {
    double origin[3] = { 0, 0, 0 };
    std::cout << "resetting the origin to center" << std::endl;
    img->SetOrigin(origin);
  }

  ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetInput(img);
  writer->SetFileName(argv[2]);
  writer->UseCompressionOn();
  writer->Write();
  return 0;
}
