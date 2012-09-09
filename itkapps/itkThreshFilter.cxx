#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#pragma warning ( disable : 4503 )
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <set>

#include <string.h>
#include <sys/types.h>
#include <stdlib.h>    // for exit, system
#include <math.h>
#include <errno.h>

#include <itkImage.h>
#include <itkImageFileReader.h> 
#include <itkImageFileWriter.h>
#include <itkExtractImageFilter.h>
#include <itkImageToImageFilter.h>
#include <itkScalarImageKmeansImageFilter.h>
#include <itkNeighborhoodIterator.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkThresholdImageFilter.h>
#include <itkCastImageFilter.h>

using namespace std;
using namespace itk;

typedef double PixelType;
typedef short  BinaryType;

enum { ImageDimension = 3 };
typedef Image<PixelType,ImageDimension>       ImageType;
typedef Image<BinaryType,ImageDimension>      BinaryImageType;
typedef ImageType::RegionType                 ImageRegionType;
typedef ImageRegionIterator< ImageType >      IteratorType;
typedef ImageRegionConstIterator<ImageType>   ConstIteratorType;
typedef NeighborhoodIterator<ImageType>       NeighborhoodIteratorType;
typedef ImageType::Pointer                    ImagePointer;

typedef ImageFileReader< ImageType >          VolumeReaderType;
typedef ImageFileWriter< ImageType >          VolumeWriterType;
typedef ImageFileWriter< BinaryImageType >    BinaryVolumeWriterType;
typedef BinaryThresholdImageFilter<ImageType,BinaryImageType> BinaryThresholdFilter;
typedef ThresholdImageFilter<ImageType> ThresholdFilter;
typedef CastImageFilter<ImageType, BinaryImageType> CastFilter;


int main(const int argc, const char **argv)
{
  if (argc < 5) {
    cout << "itkThresholdFilter input output low high [--bin]" << endl;
    return 1;
  }

  char *inputFileName = strdup(argv[1]);
  char *outputFileName = strdup(argv[2]);
  double low = atof(argv[3]);
  double high = atof(argv[4]);

	cout << "Lower threshold = " << low << endl;
	cout << "High  threshold = " << high << endl;

  int binary = 0;
  if (argc > 5 && strcmp(argv[5], "--bin") == 0) {
    binary = 1;
  }

  ImagePointer inputImage;

  VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
  imageReader->SetFileName(inputFileName) ;
  try {
    imageReader->Update();
  } catch (ExceptionObject & err) {
    cerr << "ExceptionObject caught!" << endl;
    cerr << err << endl;
    return EXIT_FAILURE;	
  }    
  inputImage = imageReader->GetOutput();
 
  BinaryVolumeWriterType::Pointer writer = BinaryVolumeWriterType::New();
  writer->SetFileName(outputFileName); 
  writer->UseCompressionOn();

  if (binary) {
    BinaryThresholdFilter::Pointer filter = BinaryThresholdFilter::New();
    filter->SetInput(inputImage);
    filter->SetLowerThreshold(low);
    filter->SetUpperThreshold(high);
    filter->SetInsideValue(1);
    filter->SetOutsideValue(0);
    filter->Update();
    writer->SetInput(filter->GetOutput());
  } else {
    ThresholdFilter::Pointer filter = ThresholdFilter::New();
    filter->SetInput(inputImage);
    filter->ThresholdOutside(low, high);
    filter->Update();
    CastFilter::Pointer caster = CastFilter::New();
    caster->SetInput(filter->GetOutput());
    caster->Update();
    writer->SetInput(caster->GetOutput());
  }
  writer->Write();
  return 0;
}
