#include <math.h>

#include <itkDefaultDynamicMeshTraits.h>
#include <itkSpatialObjectReader.h>
#include <itkSpatialObjectWriter.h>
#include <itkSceneSpatialObject.h>
#include <itkMetaMeshConverter.h>
#include <itkVTKPolyDataReader.h>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h> 
#include <itkInterpolateImageFunction.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkDanielssonDistanceMapImageFilter.h>

#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkFloatArray.h>
#include "vtkPolyDataToitkMesh.h"

#include <iostream>

typedef itk::DefaultDynamicMeshTraits<double,3,3,double,double> MeshTraitsType;
typedef itk::Mesh<double,3,MeshTraitsType> MeshType;
typedef itk::MeshSpatialObject<MeshType> MeshSOType;
typedef MeshTraitsType::PointType PointType;
typedef MeshTraitsType::CellType CellType;
typedef itk::VTKPolyDataReader<MeshType> VTKReaderType;
typedef itk::MetaMeshConverter<3,double,MeshTraitsType> MeshConverterType;
typedef itk::MeshSpatialObject<MeshType> meshSOType;

typedef itk::Image<short,3> BinaryImageType;
typedef itk::Image<double,3> ImageType;
typedef itk::InterpolateImageFunction<ImageType,double> InterpolatorType;
typedef itk::NearestNeighborInterpolateImageFunction<ImageType,double> NNInterpolatorType;
typedef itk::LinearInterpolateImageFunction<ImageType,double> LinearInterpolatorType;
typedef itk::ImageFileReader<ImageType> ImageReaderType;
typedef itk::BinaryThresholdImageFilter<ImageType,BinaryImageType> ThresholdFilterType;
typedef itk::DanielssonDistanceMapImageFilter<BinaryImageType,ImageType> DistanceFilterType;


using namespace std;

namespace {
  int extraction(char* imageName, char* meshName, char* attrName, char* refImageName) {
    cout << "Image: " << imageName << endl;
    cout << "Mesh: " << meshName << endl;
		cout << "Attribute: " << attrName << endl;
		cout << "Reference: " << refImageName << endl;
    
    DistanceFilterType::Pointer distance = DistanceFilterType::New();

    typedef itk::ImageFileReader<DistanceFilterType::VectorImageType> VectorReaderType;
    VectorReaderType::Pointer vreader = VectorReaderType::New();
    vreader->SetFileName(imageName);
    vreader->Update();
    DistanceFilterType::VectorImageType::Pointer closestVectorImage = vreader->GetOutput();

    ImageReaderType::Pointer ireader = ImageReaderType::New();
    ireader->SetFileName(refImageName);
    ireader->Update();
    ImageType::Pointer imagePtr = ireader->GetOutput();

    vtkPolyDataReader* reader = vtkPolyDataReader::New();
    reader->SetFileName(meshName);
    reader->Update();

    vtkPolyDataToitkMesh* converter = new vtkPolyDataToitkMesh();
    converter->SetInput(reader->GetOutput());
    converter->ConvertvtkToitk();
    MeshType *inputMesh = dynamic_cast<MeshType*>(converter->GetOutput());

		MeshType::PointsContainerPointer points = inputMesh->GetPoints();

    InterpolatorType::Pointer interpolator = NULL;
    interpolator = static_cast<InterpolatorType::Pointer>(NNInterpolatorType::New());
		interpolator->SetInputImage(imagePtr);
    InterpolatorType::ContinuousIndexType pointIndex;

    ImageType::SpacingType spacing = imagePtr->GetSpacing();
    ImageType::PointType origin = imagePtr->GetOrigin();

    int numPoints = points->Size();
    double* extractValue = new double[3*numPoints];
    memset(extractValue, 0, sizeof(extractValue));

    int e[3] = { 1, 1, 1 };
    ofstream attr;
    attr.open(attrName);

    for (int i = 0; i < numPoints; i++) { 
      PointType p = points->GetElement(i);
      for (unsigned int d = 0; d < 3; d++) {
        pointIndex[d] = (e[d] * p[d] - origin[d]) / spacing[d];
      }
      double ev = interpolator->EvaluateAtContinuousIndex(pointIndex);
      if (ev <= 0) {
        DistanceFilterType::VectorImageType::PixelType vectorPixel;
        DistanceFilterType::VectorImageType::IndexType vectorIndex;
        for (int d = 0; d < 3; d++) {
          vectorIndex[d] = (long int) round(pointIndex[d]);
        }
        vectorPixel = closestVectorImage->GetPixel(vectorIndex);
        for (int d = 0; d < 3; d++) {
          extractValue[3*i+d] = vectorPixel[d] * spacing[d];
        }
        attr << extractValue[3*i] << "," << extractValue[3*i+1] << "," << extractValue[3*i+2] << endl;
      } else {
        attr << "0,0,0" << endl;
      }
    }
    attr.close();

    return 0;
  }
}

int main(int argc, char* argv[]) {
  if (argc < 5) {
    cout << "usage: " << argv[0] << " distance-vector-image vtk-mesh attr-name-to-save ref-image" << endl;
    return 0;
  }

  ::extraction(argv[1], argv[2], argv[3], argv[4]);
}
