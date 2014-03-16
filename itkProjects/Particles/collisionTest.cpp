//
//  collisionTest.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 2/9/13.
//
//

#include "collisionTest.h"

#include "vtkCell.h"
#include "vtkTriangle.h"
#include "vtkPolyData.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkTriangleFilter.h"
#include "vtkPolyDataNormals.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"

#include "ozcollide/ozcollide.h"
#include "ozcollide/aabbtree_poly.h"
#include "ozcollide/aabbtreepoly_builder.h"
#include "ozcollide/vector.h"
#include "ozcollide/vec3f.h"

ozcollide::Vector<ozcollide::Vec3f> ozVerts;
ozcollide::Vector<ozcollide::Polygon> ozPolys;

int main(int argc, char* argv[]) {
    const char* file = "/data/Particles/surfaces/00.vtp";
    vtkXMLPolyDataReader* reader = vtkXMLPolyDataReader::New();
    reader->SetFileName(file);
    reader->Update();
    vtkPolyData* poly = reader->GetOutput();

    vtkTriangleFilter* triangleFilter = vtkTriangleFilter::New();
    triangleFilter->SetInput(poly);
    triangleFilter->Update();
    poly = triangleFilter->GetOutput();

    vtkPolyDataNormals* normalFilter = vtkPolyDataNormals::New();
    normalFilter->SetInput(poly);
    normalFilter->Update();
    poly = normalFilter->GetOutput();

    vtkFloatArray* normalArray = vtkFloatArray::SafeDownCast(poly->GetPointData()->GetNormals());


    const int nPoints = poly->GetNumberOfPoints();
    const int nTris = poly->GetNumberOfCells();

    vtkPoints* points = poly->GetPoints();
    for (int i = 0; i < nPoints; i++) {
        ozcollide::Vec3f vert;
        double* pts = points->GetPoint(i);
        vert.set(pts[0], pts[1], pts[2]);
        ozVerts.add(vert);
    }

    for (int i = 0; i < nTris; i++) {
        int vertIndices[3];
        vtkCell* tri = poly->GetCell(i);
        const int nVerts = tri->GetNumberOfPoints();
        for (int j = 0; j < nVerts; j++) {
            int id = tri->GetPointId(j);
            vertIndices[j] = id;
        }
        ozcollide::Vec3f normal;
        double* normalTuple = normalArray->GetTuple3(i);
        normal.set(normalTuple[0], normalTuple[1], normalTuple[2]);
        ozcollide::Polygon* p = new ozcollide::Polygon();
        p->setIndicesMemory(3, vertIndices);
        p->setNormal(normal);
        ozPolys.add(p[0]);
    }

    cout << "# points: " << poly->GetNumberOfPoints() << endl;
    cout << "# cells : " << poly->GetNumberOfCells() << endl;

    ozcollide::AABBTreePolyBuilder builder;
    ozcollide::AABBTreePoly *baseTree = builder.buildFromPolys(ozPolys.mem(), nTris, ozVerts.mem(), nPoints);

    ozPolys.clear();
    ozVerts.clear();

    reader->Delete();
    normalFilter->Delete();
    triangleFilter->Delete();
}
