#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"
#include "vtkCell.h"
#include "vtkMath.h"
#include "vtkKdTreePointLocator.h"
#include "vtkIdList.h"
#include "vtkPlane.h"
#include "vtkTriangle.h"
#include "vtkTriangleFilter.h"

#include <iostream>

using namespace std;

static void projectMeshWithLocator(vtkPolyData* m1, vtkPolyData* m2, char* filename) {
  if (m2->GetCell(0)->GetCellType() == VTK_TRIANGLE_STRIP) {
    vtkTriangleFilter* filter = vtkTriangleFilter::New();
    filter->SetInput(m2);
    filter->Update();
    m2 = filter->GetOutput();
  }

  vtkKdTreePointLocator* pointLocator = vtkKdTreePointLocator::New();
  pointLocator->SetDataSet(m2);
  pointLocator->BuildLocator();

  vtkPoints* pl1 = m1->GetPoints();
  m2->BuildLinks();

  double x[3], cp[3];
  for (int i = 0; i < pl1->GetNumberOfPoints(); i++) {
    pl1->GetPoint(i, x);
    int id = pointLocator->FindClosestPoint(x);
    m2->GetPoint(id, cp);

    vtkIdList* cellIds = vtkIdList::New();
    m2->GetPointCells(id, cellIds);
    
    /**
     * projection a point onto cell surface
     * 1) find cell
     * 2) compute normal 
     * 3) construct implicit plane function
     * 4) projection onto the plane
     * 5) determine the projected point is inside the cell
     * 6) find the point with minimum distance
     */
    double txp[3], xp[3], dist2 = vtkMath::Distance2BetweenPoints(x, cp);
    memcpy(txp, cp, sizeof(cp));

    for (int j = 0; j < cellIds->GetNumberOfIds(); j++) {
      vtkCell* cell = m2->GetCell(j);

      double t1[3], t2[3], t3[3], normal[3];
      cell->GetPoints()->GetPoint(0, t1);
      cell->GetPoints()->GetPoint(1, t2);
      cell->GetPoints()->GetPoint(2, t3);

      vtkTriangle::ComputeNormal(t1, t2, t3, normal);
      vtkPlane::GeneralizedProjectPoint(x, t1, normal, xp);

      int inside = vtkTriangle::PointInTriangle(xp, t1, t2, t3, 0);
      double ndist2 = vtkMath::Distance2BetweenPoints(x, xp);
      if (inside == 1) {
        if (dist2 > ndist2) {
          memcpy(txp, xp, sizeof(xp));
          cout << i << " vertex found closer projection onto a polygon (" << ndist2 << " < " << dist2 << ")" << endl;
          dist2 = ndist2;
        }
      } else {
        if (dist2 > ndist2) {
        }
      }
    }

   m1->GetPoints()->SetPoint(i, txp);
  }

  vtkPolyDataWriter* writer = vtkPolyDataWriter::New();
  writer->SetFileName(filename);
  writer->SetInput(m1);
  writer->Write();
}

int main(int argc, char *argv[]) {
  if (argc < 3) {
    cout << "usage: " << argv[0] << " mesh1 mesh2 mesh-out" << endl;
  }

  vtkPolyDataReader* r1 = vtkPolyDataReader::New();
  r1->SetFileName(argv[1]);
  r1->Update();

  vtkPolyDataReader* r2 = vtkPolyDataReader::New();
  r2->SetFileName(argv[2]);
  r2->Update();

  vtkPolyData* m1 = r1->GetOutput();
  vtkPolyData* m2 = r2->GetOutput();

  // projectMesh(m1, m2);
  projectMeshWithLocator(m1, m2, argv[3]);
}
