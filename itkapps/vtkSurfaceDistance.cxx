#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"

#include <iostream>

using namespace std;

// compute distance from m1 to m2
static void computeDistance(vtkPolyData* m1, vtkPolyData* m2) {
  vtkPoints* pl1 = m1->GetPoints();
  for (int i = 0; i < pl1->GetNumberOfPoints(); i++) {
    double p[3];
    pl1->GetPoint(i, p);
    cout << p[0] << "," << p[1] << "," << p[2] << endl;
  }
}

int main(int argc, char *argv[]) {
  if (argc < 3) {
    cout << "usage: " << argv[0] << " mesh1 mesh2 " << endl;
  }

  vtkPolyDataReader* r1 = vtkPolyDataReader::New();
  r1->SetFileName(argv[1]);
  r1->Update();

  vtkPolyDataReader* r2 = vtkPolyDataReader::New();
  r2->SetFileName(argv[2]);
  r2->Update();

  vtkPolyData* m1 = r1->GetOutput();
  vtkPolyData* m2 = r2->GetOutput();

  computeDistance(m1, m2);
}
