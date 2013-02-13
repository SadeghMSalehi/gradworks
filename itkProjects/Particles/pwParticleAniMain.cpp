//
//  pwParticleAniMain.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 2/11/13.
//
//

#include "pwParticleAniMain.h"

#include "vtkPolyDataReader.h"
#include "vtkPolyData.h"
#include "vtkCellTreeLocator.h"
#include "vtkTriangle.h"
#include "vtkGenericCell.h"

int main(int argc, char* argv[]) {
    vtkPolyData* mesh = NULL;
    vtkPolyDataReader* reader = vtkPolyDataReader::New();
    reader->SetFileName("/data/Particles/ellipse3use.vtk ");
    reader->Update();
    mesh = reader->GetOutput();
    mesh->Print(cout);

    vtkCellTreeLocator* loc = vtkCellTreeLocator::New();
    loc->SetDataSet(mesh);
    loc->Update();

    double p0[3] = { 40, 40, 40 };
    double p1[3] = { 0, 40, 40 };
    double x[3];
    double pCoord[3];
    double t;
    int subId;
    vtkIdType cellId;
    vtkGenericCell* cell = vtkGenericCell::New();

    loc->IntersectWithLine(p0, p1, 0, t, x, pCoord, subId, cellId, cell);
    cout << x[0] << "," << x[1] << "," << x[2] << endl;
}