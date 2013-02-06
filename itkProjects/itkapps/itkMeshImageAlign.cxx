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
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>

using namespace std;
using namespace itk;

typedef short PixelType;
typedef unsigned char BinaryType;

enum { ImageDimension = 3 };

typedef Image<PixelType,ImageDimension>       ImageType;
typedef ImageFileReader<ImageType>            VolumeReaderType;

int main(const int argc, const char **argv)
{
  //PARSE_ARGS;

  if (argc < 4) {
    cout << "usage: itkMeshImageAlign in-image in-vtk out-vtk" << endl;
    return 0;
  }

  string imageName = string(argv[1]);
  string inVtkName = string(argv[2]);
  string outVtkName = string(argv[3]);

  ImageType::Pointer inputImage ;

  // load image
  cout << "Loading file " << imageName << endl;

  VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
  imageReader->SetFileName(imageName.c_str()) ;

  try {
    imageReader->Update();
  } catch (ExceptionObject & err) {
    cerr << "ExceptionObject caught!" << endl;
    cerr << err << endl;
    return EXIT_FAILURE;	
  }    

  inputImage = imageReader->GetOutput();

  ImageType::PointType origin = inputImage->GetOrigin();

  vtkPolyDataReader* vtkReader = vtkPolyDataReader::New();
  vtkReader->SetFileName(inVtkName.c_str());
  vtkReader->Update();

  vtkTransform* transform = vtkTransform::New();
  transform->Translate(-origin[0], -origin[1], -origin[2]);

  vtkTransformPolyDataFilter* txFilter = vtkTransformPolyDataFilter::New();
  txFilter->SetInput(vtkReader->GetOutput());
  txFilter->SetTransform(transform);
  txFilter->Update();

  vtkPolyDataWriter* vtkWriter = vtkPolyDataWriter::New();
  vtkWriter->SetFileName(outVtkName.c_str());
  vtkWriter->SetInput(txFilter->GetOutput());
  vtkWriter->Write();

  return 0;
}
