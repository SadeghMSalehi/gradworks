//
//  vGraph.h
//  ktools
//
//  Created by Joowhi Lee on 8/20/15.
//
//

#ifndef __ktools__vGraph__
#define __ktools__vGraph__

#include <stdio.h>

#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkDataArray.h>

#include <unordered_set>
#include <unordered_map>
#include <vector>

struct vGraph {
    enum Color { Black = 0, Gray = 1, White = 2 };
    
    typedef std::vector<vtkIdType> vtkIdVector;
    typedef std::vector<vtkIdVector> vtkAdjVector;
    
    vtkDataSet* data;
    vtkAdjVector adjList;
    std::vector<double> dist;
    vtkIdVector pred;
    
    vGraph(vtkDataSet* d): data(d) {}

    void buildAdjacencyList();
    void computeGeodesicDistance(vtkIdVector& source);
};


#endif /* defined(__ktools__vGraph__) */
