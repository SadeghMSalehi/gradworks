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
#include <itkConnectedComponentImageFilter.h>

using namespace std;
using namespace itk;

typedef unsigned char BinaryType;
typedef short PixelType;

enum { ImageDimension = 3 };
typedef Image<BinaryType,ImageDimension>      BinaryImageType;
typedef Image<PixelType,ImageDimension>       ImageType;
typedef ImageType::Pointer                    ImagePointer;
typedef ImageType::RegionType                 ImageRegionType;
typedef ImageRegionIterator< ImageType >      IteratorType;
typedef ImageRegionConstIterator<ImageType>   ConstIteratorType;
typedef NeighborhoodIterator<ImageType>       NeighborhoodIteratorType;

typedef ImageFileReader< ImageType >          VolumeReaderType;
typedef ImageFileWriter< ImageType >          VolumeWriterType;
typedef ConnectedComponentImageFilter<ImageType,ImageType> ConnectedComponentFilterType;

int main(const int argc, const char **argv)
{
  if (argc < 3) {
    cout << argv[0] << " input output " << endl;
    return 1;
  }

  char *inputFileName = strdup(argv[1]);
  char *outputFileName = strdup(argv[2]);

  ImagePointer inputImage ;

  VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
  imageReader->SetFileName(inputFileName) ;
  try {
    imageReader->Update();
    inputImage = imageReader->GetOutput();
  } catch (ExceptionObject & err) {
    cerr << "ExceptionObject caught!" << endl;
    cerr << err << endl;
    return EXIT_FAILURE;	
  }    
 
  ConnectedComponentFilterType::Pointer filter = ConnectedComponentFilterType::New();
  filter->SetInput(inputImage);
  filter->Update();

   // when outputFileName is set
  string outFileName = "dummy";
  if (outputFileName) {
     outFileName = string(outputFileName);
  }
  
  VolumeWriterType::Pointer writer = VolumeWriterType::New();
  writer->SetFileName(outFileName.c_str()); 
  writer->SetInput(filter->GetOutput());
  writer->Write();

  return 0;
}
