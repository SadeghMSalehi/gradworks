#include "myParticleCore.h"
#include "myParticleBSpline.h"

#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridWriter.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"


#include "iostream"

using namespace std;
using namespace pi;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cout << "usage: %s particle_file.txt output.vtu" << endl;
        return 0;
    }

    ParticleSystem sys;
    sys.LoadSystem(argv[1], 0);

    const int nPoints = sys.GetNumberOfParticles();

    vtkUnstructuredGrid* grid = vtkUnstructuredGrid::New();
    vtkPoints* points = vtkPoints::New();
    points->SetNumberOfPoints(nPoints);
    for (int i = 0; i < nPoints; i++) {
        points->SetPoint(i, sys[0][i].x);
    }

    vtkDoubleArray* displacement = vtkDoubleArray::New();
    displacement->SetName("Displacement");
    displacement->SetNumberOfComponents(3);
    displacement->SetNumberOfValues(nPoints);
    for (int i = 0; i < nPoints; i++) {
        displacement->SetTuple3(i,
                                sys[1][i].x[0] - sys[0][i].x[0],
                                sys[1][i].x[1] - sys[0][i].x[1],
                                sys[1][i].x[2] - sys[0][i].x[2]);
    }

    grid->SetPoints(points);
    grid->GetPointData()->SetVectors(displacement);

    vtkUnstructuredGridWriter* writer = vtkUnstructuredGridWriter::New();
    writer->SetFileName(argv[2]);
#if VTK_MAJOR_VERSION >= 6
    writer->SetInputData(grid);
#else
    writer->SetInput(grid);
#endif
    writer->Write();
    
    return 0;
}