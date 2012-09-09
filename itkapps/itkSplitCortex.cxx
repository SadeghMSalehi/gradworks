#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#pragma warning ( disable : 4503 )
#endif

#include <string>
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
#include <itkArray.h>
#include <itkImageRegionIteratorWithIndex.h>
#include "itkPlaneSpatialFunction.h"
#include "itkSplitCortexCLP.h"

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
typedef BinaryImageType::Pointer              BinaryImagePointer;

typedef ImageFileReader< ImageType >          VolumeReaderType;
typedef ImageFileWriter< ImageType >          VolumeWriterType;

typedef ImageFileReader< BinaryImageType >    BinaryVolumeReaderType;
typedef ImageFileWriter< BinaryImageType >    BinaryVolumeWriterType;
typedef ExtractImageFilter<BinaryImageType,BinaryImageType> ExtractImageFilterType;

static int debug = 1;

int main(const int argc, const char **argv)
{
  //PARSE_ARGS;

  if (argc < 9) {
    cout << "usage: itkSplitCortex in out no nx ny nz ox oy oz" << endl;
    return 0;
  }

  string imageName = string(argv[1]);
  string outImageName = string(argv[2]);
  int planeNo = atoi(argv[3]);

  double n[3] = { 1, 0, -0.002 };
  if (argc >= 6) {
    for (int i = 0; i < 3; i++) {
      n[i] = atof(argv[i+4]);
    }
  }


  BinaryImagePointer inputImage ;

  // load image
  if (debug) {
    cout << "Loading file " << imageName << endl;
  }

  BinaryVolumeReaderType::Pointer imageReader = BinaryVolumeReaderType::New();
  imageReader->SetFileName(imageName.c_str()) ;
  try {
    imageReader->Update();
  } catch (ExceptionObject & err) {
    cerr << "ExceptionObject caught!" << endl;
    cerr << err << endl;
    return EXIT_FAILURE;	
  }    
  inputImage = imageReader->GetOutput();

  typedef itk::ImageRegionIteratorWithIndex<BinaryImageType> BinaryIteratorType;
  BinaryIteratorType it(inputImage, inputImage->GetRequestedRegion());
  BinaryIteratorType::IndexType imageIndex;

  typedef itk::PlaneSpatialFunction<> PlaneFunctionType;
  PlaneFunctionType::Pointer planeFunc = PlaneFunctionType::New();
  BinaryImageType::SizeType inputImageSize;
  inputImageSize = inputImage->GetLargestPossibleRegion().GetSize();

  Point<double,3> origin;
  for (int i = 0; i < 3; i++) {
    origin[i] = inputImageSize[i] / 2.0;
  }

  if (argc >= 10) {
    for (int i = 0; i < 3; i++) {
      origin[i] = atoi(argv[i+7]);
    }
  }

  planeFunc->SetOrigin(origin);
  planeFunc->SetNormal(n[0], n[1], n[2]);

  cout << "Normal: " << n[0] << "," << n[1] << "," << n[2] << endl;
  cout << "Origin: " << origin[0] << "," << origin[1] << "," << origin[2] << endl;

  Point<double,3> pos;

  int voxelCounter = 0;
  int maxVoxels = inputImageSize[0] * inputImageSize[1] * inputImageSize[2];
  int percent = 0;

  for (;!it.IsAtEnd();++it) {
    voxelCounter++;

    imageIndex = it.GetIndex();
    pos[0] = imageIndex[0];
    pos[1] = imageIndex[1];
    pos[2] = imageIndex[2];

    if (it.Get() != 3 && it.Get() != 4) {
      continue;
    }

    bool inOut = planeFunc->Evaluate(pos);

    if (inOut) {
      it.Set(3);
    } else {
      it.Set(4);
    }
  }

  cout << "Writing " << outImageName << endl;
  BinaryVolumeWriterType::Pointer writer = BinaryVolumeWriterType::New();
  writer->SetFileName(outImageName.c_str()); 
  writer->SetInput(inputImage);
  writer->UseCompressionOn();
  writer->Write();
  return 0;
}
