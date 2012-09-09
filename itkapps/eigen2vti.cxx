#include "iostream"
#include "itkImage.h"
#include "itkImageCommon.h"
#include "itkImageRegion.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "vtkUniformGrid.h"
#include "vtkXMLImageDataWriter.h"
#include "vtkPointData.h"
#include "vtkDoubleArray.h"

typedef itk::Vector<double,3> VectorType;
typedef itk::Image<VectorType,3> VectorImageType;
typedef itk::Image<double,3> ScalarImageType;


int main(int argc, char* argv[]) {
  if (argc < 4) {
    cout << "usage: " << argv[0] << " input-image vector-image output-vti" << endl;
    return 0;
  }
  
  int ret = 0;
  ScalarImageType::Pointer sPtr = ReadImageT<ScalarImageType>(argv[1], ret);
  ScalarImageType::RegionType region = sPtr->GetRequestedRegion();
  ScalarImageType::SizeType size = region.GetSize();
  ScalarImageType::SpacingType spacing = sPtr->GetSpacing();
  ScalarImageType::PointType origin = sPtr->GetOrigin();
  
  VectorImageType::Pointer vPtr = ReadImageT<VectorImageType>(argv[2], ret);

  vtkUniformGrid* gridPtr = vtkUniformGrid::New();
  gridPtr->SetDimensions(size[0], size[1], size[2]);
  gridPtr->SetScalarTypeToDouble();
  gridPtr->SetSpacing(spacing[0], spacing[1], spacing[2]);
  gridPtr->SetOrigin(origin[0], origin[1], origin[2]);
  gridPtr->AllocateScalars();

  cout << "Number of voxels: " << gridPtr->GetNumberOfPoints() << endl;

  vtkDoubleArray* eigenVector = vtkDoubleArray::New();
  eigenVector->SetNumberOfComponents(3);
  eigenVector->SetNumberOfTuples(gridPtr->GetNumberOfPoints());
  eigenVector->SetName("Eigen Vector");

  gridPtr->GetPointData()->SetVectors(eigenVector);

  typedef itk::ImageRegionIteratorWithIndex<ScalarImageType> ScalarIteratorType;
  typedef itk::ImageRegionIteratorWithIndex<VectorImageType> VectorIteratorType;

  ScalarIteratorType sIter(sPtr, region);
  VectorIteratorType vIter(vPtr, region);

  double* scalars = (double*) gridPtr->GetScalarPointer();
  int i = 0;
  for (; !sIter.IsAtEnd() && !vIter.IsAtEnd(); ++sIter, ++vIter, ++i) {
    double s = sIter.Get();
    VectorType v = vIter.Get();
    double vd[3];
    vd[0] = v[0]; vd[1] = v[1]; vd[2] = v[2];

    scalars[i] = s;
    eigenVector->SetTuple(i, vd);
  }

  cout << "Writing image file: " << argv[3] << endl;
  vtkXMLImageDataWriter* vtiWriter = vtkXMLImageDataWriter::New();
  vtiWriter->SetInput(gridPtr);
  vtiWriter->SetFileName(argv[3]);
  vtiWriter->Write();
}
