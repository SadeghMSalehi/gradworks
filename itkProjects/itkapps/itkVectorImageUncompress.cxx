#include <iostream>

#include <itkVector.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

typedef itk::Vector<float,3> VectorType;
typedef itk::Image<VectorType,3> ImageType;
typedef itk::ImageFileReader<ImageType> ImageReaderType;
typedef itk::ImageFileWriter<ImageType> ImageWriterType;

int main(int argc, char* argv[]) {
  if (argc < 3) {
    std::cout << "usage: " << argv[0] << " image-in image-out" << std::endl;
  }

  ImageReaderType::Pointer reader = ImageReaderType::New();
  reader->SetFileName(argv[1]);
  reader->Update();

  ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetFileName(argv[2]);
  writer->SetInput(reader->GetOutput());
  writer->UseCompressionOff();
  writer->Write();
}
