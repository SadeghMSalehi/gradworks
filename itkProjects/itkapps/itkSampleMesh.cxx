#include "itkSampleMeshCLP.h"

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
  int extraction(string image, string mesh, string intp, string dstv, double lowThreshold, string attrfileName) {
    cout << "Image: " << image << endl;
    cout << "Mesh: " << mesh << endl;
		cout << "Interpolation: " << intp << endl;
    
    ImageReaderType::Pointer ireader = ImageReaderType::New();
    ireader->SetFileName(image.c_str());
    ireader->Update();
    ImageType::Pointer imagePtr = ireader->GetOutput();

    InterpolatorType::Pointer interpolator = NULL;
		if (intp == "nn") {
			interpolator = static_cast<InterpolatorType::Pointer>(NNInterpolatorType::New());
		} else if (intp == "linear") {
			interpolator = static_cast<InterpolatorType::Pointer>(LinearInterpolatorType::New());
		}
		interpolator->SetInputImage(imagePtr);
    InterpolatorType::ContinuousIndexType pointIndex;

    /*
		MeshConverterType* converter = new MeshConverterType();
		MeshSOType::Pointer soMesh = converter->ReadMeta(mesh.c_str());
		MeshType::Pointer inputMesh = soMesh->GetMesh();
		MeshType::PointsContainerPointer points = inputMesh->GetPoints();
    */
    vtkPolyDataReader* reader = vtkPolyDataReader::New();
    reader->SetFileName(mesh.c_str());
    reader->Update();
    //vtkPolyData* vtkMesh = reader->GetOutput();

    vtkPolyDataToitkMesh converter;
    converter.SetInput(reader->GetOutput());
    converter.ConvertvtkToitk();
    MeshType *inputMesh = dynamic_cast<MeshType*>(converter.GetOutput());
		MeshType::PointsContainerPointer points = inputMesh->GetPoints();

    DistanceFilterType::Pointer distance = DistanceFilterType::New();
    DistanceFilterType::VectorImageType::Pointer closestVectorImage = NULL;

    bool dstvExists = itksys::SystemTools::FileExists(dstv.c_str(), true);

    if (!dstvExists) {
      cout << "Generating distance vector: " << dstv << endl;
      ThresholdFilterType::Pointer threshold = ThresholdFilterType::New();

      threshold->SetInput(imagePtr);
      threshold->SetLowerThreshold(1);
      threshold->SetUpperThreshold(32000);

      distance->UseImageSpacingOn();
      distance->SquaredDistanceOn();
      distance->SetInput(threshold->GetOutput());
      distance->Update();

      closestVectorImage = distance->GetVectorDistanceMap();

      typedef itk::ImageFileWriter<DistanceFilterType::VectorImageType> VectorWriterType;

      VectorWriterType::Pointer writer = VectorWriterType::New();
      writer->SetInput(closestVectorImage);
      writer->UseCompressionOn();
      writer->SetFileName(dstv);
      writer->Write();

      cout << "Generating distance vector done" << endl;
    } else {
      cout << "Loading distance vector: " << dstv << endl;
      typedef itk::ImageFileReader<DistanceFilterType::VectorImageType> VectorReaderType;

      VectorReaderType::Pointer vreader = VectorReaderType::New();
      vreader->SetFileName(dstv);
      vreader->Update();

      closestVectorImage = vreader->GetOutput();
    }

    ImageType::SpacingType spacing = imagePtr->GetSpacing();
    ImageType::PointType origin = imagePtr->GetOrigin();

    int numPoints = points->Size();
    double* extractValue = new double[numPoints];
    memset(extractValue, 0, sizeof(extractValue));

    int e[3] = { 1, 1, 1 };
    ofstream attr;
    attr.open(attrfileName.c_str());

    cout << "Sampling points... " << endl;
    for (int i = 0; i < numPoints; i++) { 
      PointType p = points->GetElement(i);
      for (unsigned int d = 0; d < 3; d++) {
        pointIndex[d] = (e[d] * p[d] - origin[d]) / spacing[d];
      }
      extractValue[i] = interpolator->EvaluateAtContinuousIndex(pointIndex);
      // cout << "(" << pointIndex[0] << ", " << pointIndex[1] << ", " << pointIndex[2] << " ) => (";
      if (extractValue[i] <= lowThreshold) {
        DistanceFilterType::VectorImageType::PixelType vectorPixel;
        DistanceFilterType::VectorImageType::IndexType vectorIndex;
        for (int d = 0; d < 3; d++) {
          vectorIndex[d] = (long int) round(pointIndex[d]);
        }
        vectorPixel = closestVectorImage->GetPixel(vectorIndex);
        for (int d = 0; d < 3; d++) {
          pointIndex[d] = vectorIndex[d] + vectorPixel[d];
        }
        extractValue[i] = interpolator->EvaluateAtContinuousIndex(pointIndex); 
        if (extractValue[i] <= lowThreshold) {
          cerr << "error at extracting closest point " << endl;
        }
      }
      // cout << pointIndex[0] << "," << pointIndex[1] << "," << pointIndex[2] << ") " << endl;
      attr << extractValue[i] << endl;
    }
    attr.close();

    cout << "sampling points ... done" << endl;

    return 0;
  }
}

int main(int argc, char* argv[]) {
  PARSE_ARGS;

  itksys::SystemTools::ChangeDirectory(workingDirectory.c_str());
  ::extraction(imageName, meshName, interpolation, distanceVector, lowThreshold, attrfileName);
}
