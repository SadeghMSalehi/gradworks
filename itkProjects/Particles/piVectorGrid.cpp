//
//  piVectorGrid.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 2/13/13.
//
//

#include "piVectorGrid.h"

#include "vtkUnstructuredGrid.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkDoubleArray.h"
#include "vtkUnstructuredGridWriter.h"

namespace pi {
    VectorGrid::VectorGrid() {
        m_Grid = vtkUnstructuredGrid::New();
        m_Points = vtkPoints::New();
    }

    VectorGrid::~VectorGrid() {
        m_Grid->Delete();
    }

    int VectorGrid::AddPoint(double *x) {
        return m_Points->InsertNextPoint(x);
    }


    int VectorGrid::CreateScalar(const char *name) {
        m_Array = vtkDoubleArray::New();
        m_Array->SetName(name);
        m_Array->SetNumberOfComponents(1);
        return m_Grid->GetPointData()->AddArray(m_Array);
    }

    int VectorGrid::CreateVector(const char *name) {
        m_Array = vtkDoubleArray::New();
        m_Array->SetName(name);
        m_Array->SetNumberOfComponents(3);
        return m_Grid->GetPointData()->AddArray(m_Array);
    }

    int VectorGrid::AddScalar(const char *name, double s) {
        vtkDoubleArray* arr = vtkDoubleArray::SafeDownCast(m_Grid->GetPointData()->GetArray(name));
        return arr->InsertNextValue(s);

    }
    
    int VectorGrid::AddVector(const char *name, double *v) {
        vtkDoubleArray* arr = vtkDoubleArray::SafeDownCast(m_Grid->GetPointData()->GetArray(name));
        return arr->InsertNextTupleValue(v);
    }

    void VectorGrid::WriteFile(const char *filename) {
        m_Grid->SetPoints(m_Points);

        vtkUnstructuredGridWriter* w = vtkUnstructuredGridWriter::New();
        w->SetFileName(filename);
        w->SetInput(m_Grid);
        w->Update();
        w->Delete();
    }

}