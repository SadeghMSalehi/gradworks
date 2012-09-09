#include <iostream>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include "itkVectorConnectedImageFilter.h"

using namespace std;

typedef itk::Image<unsigned short,3> ShortImageType;
typedef itk::Image<unsigned int, 3> IntImageType;
//typedef itk::Vector<float,3> VectorType;
typedef itk::Image<VectorType,3> VectorImageType;
typedef itk::ImageFileReader<ShortImageType> ShortImageReaderType;
typedef itk::ImageFileReader<VectorImageType> VectorImageReaderType;
typedef itk::ImageFileWriter<ShortImageType> ShortImageWriterType;
typedef itk::ImageFileWriter<IntImageType> IntImageWriterType;
typedef itk::ImageRegionConstIteratorWithIndex<ShortImageType> SeedIteratorType;
typedef itk::VectorConnectedImageFilter<VectorImageType,IntImageType> VectorSegFilter;

int main(int argc, char* argv[]) {
  if (argc < 6) {
    cout << "usage: " << argv[0] << " dti-file seed-file output-file l1, u1, l2, u2, ..." << endl;
    exit(0);
  }

  ShortImageReaderType::Pointer seedReader = ShortImageReaderType::New();
  seedReader->SetFileName(argv[2]);
  seedReader->Update();
  ShortImageType::Pointer seedImage = seedReader->GetOutput();

  VectorImageReaderType::Pointer vectorReader = VectorImageReaderType::New();
  vectorReader->SetFileName(argv[1]);
  vectorReader->Update();
  VectorImageType::Pointer vectorImage = vectorReader->GetOutput();

  SeedIteratorType it(seedImage, seedImage->GetRequestedRegion());
  it.GoToBegin();

  VectorSegFilter::Pointer filter = VectorSegFilter::New();
  filter->SetInput(vectorImage);

  unsigned short maxLabel = 1;
  filter->AddLower(atof(argv[4]));
  filter->AddUpper(atof(argv[5]));
	filter->AddRatio(1.0);

  for (it = it.Begin(); !it.IsAtEnd(); ++ it) {
    ShortImageType::IndexType idx = it.GetIndex();
    if (it.Get() == 1) {
      filter->AddSeed(idx);
      cout << "1) " << idx << endl;
    } else {
      if (maxLabel < it.Get()) {
        maxLabel = it.Get();
      }
    }
  }

  cout << "Max Label: " << maxLabel << endl;

  for (int i = 2; i <= maxLabel; i++) {
    it.GoToBegin();
    filter->BeginNextSeed();

    filter->AddLower(atof(argv[4 + (i-1)*2]));
    filter->AddUpper(atof(argv[5 + (i-1)*2]));
    filter->AddRatio(1.0);

    for (it = it.Begin(); !it.IsAtEnd(); ++ it) {
      ShortImageType::IndexType idx = it.GetIndex();
      if (it.Get() == i) {
          filter->AddSeed(idx);
          cout << i << ") " << idx << endl;
      }
    }
  }


  filter->Update();
  IntImageType::Pointer vectorSeg = filter->GetOutput();
  IntImageWriterType::Pointer writer = IntImageWriterType::New();
  writer->SetFileName(argv[3]);
  writer->SetInput(vectorSeg);
  writer->UseCompressionOn();
  writer->Update();

}
