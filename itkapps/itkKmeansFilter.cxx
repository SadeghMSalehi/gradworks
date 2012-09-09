#include <iostream>

#include <itkImage.h>
#include <itkScalarImageKmeansImageFilter.h>
#include <itkArray.h>
#include "itkImageCommon.h"

using namespace std;

typedef itk::Image<double,3> ImageType;
typedef itk::Image<unsigned char,3>  BinaryImageType;
typedef itk::ScalarImageKmeansImageFilter<ImageType> KmeansFilterType;
typedef itk::Array<double> ClassMeanArray;

int main(const int argc, const char **argv)
{
  if (argc < 4) {
    cout << "itkKmeansFilter input output [class1 mean] [class2 mean] [class3 mean] ... " << endl;
    return 1;
  }

	int ret = 0;
  ImageType::Pointer inputImage = ReadImageT<ImageType>(argv[1], ret);

  KmeansFilterType::Pointer kmeansFilter = KmeansFilterType::New();
  kmeansFilter->SetInput(inputImage);
  kmeansFilter->UseNonContiguousLabelsOff();

  for (int i = 3; i < argc; i++) {
    cout << "Class: " << (i - 2) << ": " << static_cast<ImageType::PixelType>(atof(argv[i])) << endl;
    kmeansFilter->AddClassWithInitialMean(static_cast<ImageType::PixelType>(atof(argv[i])));
  }
  try {
    kmeansFilter->Update();
  } catch (itk::ExceptionObject& err) {
    cerr << "ExceptionObject caught!" << endl;
    cerr << err << endl;
    return EXIT_FAILURE;	
  }    

  ClassMeanArray finalMeans = kmeansFilter->GetFinalMeans();
  for (unsigned int i = 0; i < finalMeans.GetSize(); i++) {
    cout << "Class: " << (i + 1) << " = " << finalMeans[i] << endl;
  }

	WriteImageT<BinaryImageType>(argv[2], kmeansFilter->GetOutput());
  return 0;
}
