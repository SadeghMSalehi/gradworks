#include "itkImageCommon.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkDoubleArray.h"

#include "vtkStructuredGrid.h"
#include "vtkXMLStructuredGridWriter.h"

#include "iostream"

using namespace std;

typedef itk::Vector<double,3> VectorType;
typedef itk::Image<VectorType,3> VectorImageType;
typedef itk::ImageRegionIteratorWithIndex<VectorImageType> IteratorType;

int main(int argc, char* argv[]) {
  if (argc < 3) {
    cout << "usage: " << argv[0] << " input-image output-vti" << endl;
    return 0;
  }

  int ret = 0;
  VectorImageType::Pointer srcImg = ReadImageT<VectorImageType>(argv[1], ret);
  VectorImageType::RegionType srcRegion = srcImg->GetRequestedRegion();
  VectorImageType::SizeType srcSize = srcRegion.GetSize();
  
  vtkStructuredGrid* vtkGrid = vtkStructuredGrid::New();
  vtkGrid->SetDimensions(srcSize[0], srcSize[1], srcSize[2]);
  vtkGrid->SetExtent(0, srcSize[0] - 1, 0, srcSize[1] - 1, 0, srcSize[2] - 1);
  vtkGrid->Initialize();

  vtkDoubleArray* eigenVectors = vtkDoubleArray::New();
  eigenVectors->SetNumberOfComponents(3);
  eigenVectors->SetName("EigenVectors");
  
  vtkPoints* points = vtkPoints::New();
  points->SetNumberOfPoints(srcSize[0] * srcSize[1] * srcSize[2]);

  IteratorType iter(srcImg, srcRegion);
  for (int i = 0; !iter.IsAtEnd(); ++iter, i++) {
    VectorImageType::IndexType idx = iter.GetIndex();
    points->SetPoint(0, idx[0], idx[1], idx[2]);

    VectorType v = iter.Get();
    eigenVectors->InsertNextTuple3(v[0], v[1], v[2]);
  }
  vtkGrid->SetPoints(points);
  vtkGrid->GetPointData()->SetVectors(eigenVectors);

  vtkXMLStructuredGridWriter* gridWriter = vtkXMLStructuredGridWriter::New();
  gridWriter->SetInput(vtkGrid);
  gridWriter->SetFileName(argv[2]);
  gridWriter->SetDataModeToAppended();
  gridWriter->SetCompressorTypeToZLib();
  gridWriter->EncodeAppendedDataOff();
  gridWriter->Write();

  return 0;
}
