#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>

#include "vtkPolyDataConnectivityFilter2.h"

#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {
  if (argc < 5) {
    cout << "usage: " << argv[0] << " input-vtk output-vtk scalar-name min max" << endl;
    return 0;
  }

  vtkPolyDataReader* reader = vtkPolyDataReader::New();
  reader->SetFileName(argv[1]);
  reader->Update();
  
  vtkPolyData* inmesh = reader->GetOutput();
  inmesh->GetPointData()->SetActiveScalars(argv[3]);
  
  vtkPolyDataConnectivityFilter2* filter = vtkPolyDataConnectivityFilter2::New();
  filter->SetInput(inmesh);
  filter->ScalarConnectivityOn();
  filter->ExclusiveScalarConnectivityOn();
  filter->SetScalarRange(atof(argv[4]), atof(argv[5]));
  filter->SetExtractionModeToLargestRegion();
  filter->Update();

  vtkPolyData* outmesh = filter->GetOutput();

  vtkPolyDataWriter* writer = vtkPolyDataWriter::New();
  writer->SetFileName(argv[2]);
  writer->SetInput(outmesh);
  writer->Write();
}
