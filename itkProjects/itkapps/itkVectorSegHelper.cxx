#include <iostream>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include "itkVectorConnectedImageFilter.h"
#include "itkMinMaxVectorProductImageFilter.h"

using namespace std;

typedef itk::Image<short,3> ShortImageType;
//typedef itk::Vector<float,3> VectorType;
typedef itk::Image<VectorType,3> VectorImageType;
typedef itk::ImageFileReader<ShortImageType> ShortImageReaderType;
typedef itk::ImageFileReader<VectorImageType> VectorImageReaderType;
typedef itk::ImageFileWriter<ShortImageType> ShortImageWriterType;
typedef itk::ImageRegionConstIteratorWithIndex<ShortImageType> SeedIteratorType;
typedef itk::VectorConnectedImageFilter<VectorImageType,ShortImageType> VectorSegFilter;
typedef itk::MinMaxVectorProductImageFilter<VectorImageType,ShortImageType> MinMaxFilter;

int main(int argc, char* argv[]) {
  if (argc < 3) {
    cout << "usage: " << argv[0] << " vector-file out-file" << endl;
    exit(0);
  }

/*
  ShortImageReaderType::Pointer seedReader = ShortImageReaderType::New();
  seedReader->SetFileName(argv[2]);
  seedReader->Update();
  ShortImageType::Pointer seedImage = seedReader->GetOutput();
*/
  VectorImageReaderType::Pointer vectorReader = VectorImageReaderType::New();
  vectorReader->SetFileName(argv[1]);
  vectorReader->Update();
  VectorImageType::Pointer vectorImage = vectorReader->GetOutput();

/*
  SeedIteratorType it(seedImage, seedImage->GetRequestedRegion());
  it.GoToBegin();

  VectorSegFilter::Pointer filter = VectorSegFilter::New();
  filter->SetInput(vectorImage);

  for (it = it.Begin(); !it.IsAtEnd(); ++ it) {
    ShortImageType::IndexType idx = it.GetIndex();
    if (it.Get() != 0) {
      cout << idx << endl;
      filter->AddSeed(idx);
    }
  }

  filter->SetLower(0.5);
  filter->SetUpper(1);

  filter->Update();
*/

  MinMaxFilter::Pointer filter = MinMaxFilter::New();
	filter->SetInput(vectorImage);
	filter->Update();


  ShortImageType::Pointer vectorSeg = filter->GetOutput();
  ShortImageWriterType::Pointer writer = ShortImageWriterType::New();
  writer->SetFileName(argv[2]);
  writer->SetInput(vectorSeg);
  writer->Update();

}
