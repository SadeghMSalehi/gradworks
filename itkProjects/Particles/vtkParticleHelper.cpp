#include "vtkParticleHelper.h"

#include "vtkPoints.h"
#include "vtkLoopSubdivisionFilter.h"
#include "vtkDelaunay2D.h"
#include "vtkDelaunay3D.h"
#include "vtkUnstructuredGrid.h"
#include "vtkTriangleFilter.h"

vtkParticleHelper::vtkParticleHelper()
{
}


bool vtkParticleHelper::main(pi::Options &opts, pi::StringVector &args) {
    using namespace std;
    string testName = opts.GetString("--vtk", "");
    if (testName == "split") {
        testSplitParticles();
        return true;
    }
    return false;
}


void vtkParticleHelper::testSplitParticles() {
    pi::ParticleSubject subj(0, 4);

    subj[0].x[0] = 0; subj[0].x[1] = 0;
    subj[1].x[0] = 0; subj[1].x[1] = 1;
    subj[2].x[0] = 1; subj[2].x[1] = 1;
    subj[3].x[0] = 1; subj[3].x[1] = 0;

    splitParticles(subj);

    for (int i = 0; i < subj.GetNumberOfPoints(); i++) {
        cout << "p[" << i << "]: " << subj[i].x[0] << "," << subj[i].x[1] << endl;
    }
}


void vtkParticleHelper::splitParticles(pi::ParticleSubject& subj) {
    vtkPoints* points = vtkPoints::New();
    points->SetNumberOfPoints(subj.GetNumberOfPoints());

    vtkPolyData* inputPoly = vtkPolyData::New();
    inputPoly->SetPoints(points);

    vtkDataSet* grid = NULL;
    if (DIMENSIONS == 3) {
        for (int i = 0; i < subj.GetNumberOfPoints(); i++) {
            points->SetPoint(i, subj[i].x[0], subj[i].x[1], subj[i].x[2]);
        }
        vtkDelaunay3D* filter = vtkDelaunay3D::New();
        filter->SetInput(inputPoly);
        filter->Update();
        grid = filter->GetOutput();
    } else {
        for (int i = 0; i < subj.GetNumberOfPoints(); i++) {
            points->SetPoint(i, subj[i].x[0], subj[i].x[1], 0);
        }
        vtkDelaunay2D* filter = vtkDelaunay2D::New();
        filter->SetInput(inputPoly);
        filter->Update();
        grid = filter->GetOutput();
    }

    vtkTriangleFilter* triangle = vtkTriangleFilter::New();
    triangle->SetInput(grid);

    vtkLoopSubdivisionFilter* subdivFilter = vtkLoopSubdivisionFilter::New();
    subdivFilter->SetInputConnection(triangle->GetOutputPort());
    subdivFilter->Update();
    vtkPolyData* output = subdivFilter->GetOutput();

    const int newPoints = output->GetNumberOfPoints();
    subj.m_Particles.resize(newPoints);

    for (int i = 0; i < output->GetNumberOfPoints(); i++) {
        double* p = output->GetPoint(i);
        subj[i].Zero();
        fordim (k) {
            subj[i].x[k] = p[k];
        }
    }


    inputPoly->Delete();
    triangle->Delete();
    subdivFilter->Delete();
    points->Delete();
}
