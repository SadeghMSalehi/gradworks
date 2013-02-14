//
//  piVectorGrid.h
//  ParticlesGUI
//
//  Created by Joohwi Lee on 2/13/13.
//
//

#ifndef __ParticlesGUI__piVectorGrid__
#define __ParticlesGUI__piVectorGrid__

#include <iostream>
class vtkUnstructuredGrid;
class vtkPoints;
class vtkDoubleArray;

namespace pi {
    class VectorGrid {
    public:
        VectorGrid();
        ~VectorGrid();
        int AddPoint(double *x);

        int CreateScalar(const char* name);
        int CreateVector(const char* name);

        int AddVector(const char* name, double* v);
        int AddScalar(const char* name, double s);

        void WriteFile(const char* filename);
    private:
        vtkUnstructuredGrid* m_Grid;
        vtkPoints* m_Points;
        vtkDoubleArray* m_Array;
    };
}
#endif /* defined(__ParticlesGUI__piVectorGrid__) */
