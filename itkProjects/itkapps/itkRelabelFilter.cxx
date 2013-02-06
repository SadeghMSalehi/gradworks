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

using namespace std;
using namespace itk;

typedef short PixelType;
typedef unsigned char BinaryType;

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
typedef ScalarImageKmeansImageFilter<ImageType> KmeansFilterType;
typedef RelabelComponentImageFilter<ImageType, ImageType> RelabelFilterType;

static int debug = 1;

int main(const int argc, const char **argv)
{
  if (argc < 3) {
    cout << "itkRelabelsFilter input output " << endl;
    return 1;
  }

  char *inputFileName = strdup(argv[1]);
  char *outputFileName = strdup(argv[2]);

  ImagePointer inputImage ;

  // load image
  if (debug) {
    cout << "Loading file " << inputFileName << endl;
  }

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
 
  RelabelFilterType::Pointer relabelFilter = RelabelFilterType::New();
  relabelFilter->SetInput(inputImage);
  relabelFilter->Update();

  // write image
  if (debug) {
    cout << "writing output data " << string(outputFileName) << endl;
  }

   // when outputFileName is set
  string outFileName = "dummy";
  if (outputFileName) {
     outFileName = string(outputFileName);
  }
  
  VolumeWriterType::Pointer writer = VolumeWriterType::New();
  writer->SetFileName(outFileName.c_str()); 
  writer->SetInput(relabelFilter->GetOutput());
  writer->Write();
  return 0;
}
