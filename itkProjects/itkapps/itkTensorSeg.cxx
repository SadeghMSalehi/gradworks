#include <iostream>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkDiffusionTensor3D.h>
#include <itkImageRegionConstIteratorWithIndex.h>
#include "itkTensorConnectedImageFilter.h"

using namespace std;

typedef itk::Image<short,3> ShortImageType;
typedef itk::DiffusionTensor3D<double> TensorType;
typedef itk::Image<TensorType,3> TensorImageType;
typedef itk::ImageFileReader<ShortImageType> ShortImageReaderType;
typedef itk::ImageFileReader<TensorImageType> TensorImageReaderType;
typedef itk::ImageFileWriter<ShortImageType> ShortImageWriterType;
typedef itk::ImageRegionConstIteratorWithIndex<ShortImageType> SeedIteratorType;
typedef itk::TensorConnectedImageFilter<TensorImageType,ShortImageType> TensorSegFilter;

int main(int argc, char* argv[]) {
  if (argc < 4) {
    cout << "usage: " << argv[0] << " dti-file seed-file output-file" << endl;
    exit(0);
  }

  ShortImageReaderType::Pointer seedReader = ShortImageReaderType::New();
  seedReader->SetFileName(argv[2]);
  seedReader->Update();
  ShortImageType::Pointer seedImage = seedReader->GetOutput();

  TensorImageReaderType::Pointer tensorReader = TensorImageReaderType::New();
  tensorReader->SetFileName(argv[1]);
  tensorReader->Update();
  TensorImageType::Pointer tensorImage = tensorReader->GetOutput();

  SeedIteratorType it(seedImage, seedImage->GetRequestedRegion());
  it.GoToBegin();

  TensorSegFilter::Pointer filter = TensorSegFilter::New();
  filter->SetInput(tensorImage);

  for (it = it.Begin(); !it.IsAtEnd(); ++ it) {
    ShortImageType::IndexType idx = it.GetIndex();
    if (it.Get() != 0) {
      cout << idx << endl;
      filter->AddSeed(idx);
    }
  }

  filter->SetLower(0.85);
  filter->SetUpper(1);

  cout << filter << endl;

  filter->Update();
  ShortImageType::Pointer tensorSeg = filter->GetOutput();
  ShortImageWriterType::Pointer writer = ShortImageWriterType::New();
  writer->SetFileName(argv[3]);
  writer->SetInput(tensorSeg);
  writer->Update();

}
