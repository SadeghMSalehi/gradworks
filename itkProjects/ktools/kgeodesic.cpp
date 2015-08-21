//
//  kgeodesic.cpp
//  ktools
//
//  Created by Joowhi Lee on 8/20/15.
//
//

#include "kgeodesic.h"
#include "vtkio.h"

#include <vtkDoubleArray.h>
#include <vtkUnstructuredGrid.h>

#include <vector>
#include <deque>

using namespace pi;
using namespace std;

static vtkIO vio;


vtkPolygon* computeLocalTangentMap(vtkPolyData* g, vtkIdType p) {
    vtkNew<vtkIdList> cellIds, ptIds;
    g->GetPointCells(p, cellIds.GetPointer());
    
    vector<deque<vtkIdType> > neighbors(cellIds->GetNumberOfIds());
    
    for (size_t j = 0; j < cellIds->GetNumberOfIds(); j++) {
        vtkIdType c = cellIds->GetId(j);
        g->GetCellPoints(c, ptIds.GetPointer());
        for (size_t k = 0; k < ptIds->GetNumberOfIds(); k++) {
            vtkIdType q = ptIds->GetId(k);
            neighbors[j].push_back(q);
        }
        while (neighbors[j].front() != p) {
            vtkIdType n = neighbors[j].front();
            neighbors[j].pop_front();
            neighbors[j].push_back(n);
        }
    }
    
//    for (size_t j = 0; j < neighbors.size(); j++) {
//        for (size_t k = 0; k < neighbors[j].size(); k++) {
//            cout << neighbors[j][k] << ",";
//        }
//        cout << endl;
//    }
    
    vtkPolygon* tri = vtkPolygon::New();
    vtkIdList* triIds = tri->GetPointIds();
    vtkIdType firstId = neighbors[0][1];
    vtkIdType lastId = -1;
    for (size_t j = 1; j < neighbors[0].size(); j++) {
        lastId = neighbors[0][j];
        triIds->InsertNextId(lastId);
    }
    
    for (size_t j = 0; j < neighbors.size() && lastId != firstId; j++) {
        for (size_t k = 0; k < neighbors.size(); k++) {
            if (neighbors[k][1] == lastId) {
                lastId = neighbors[k][2];
                if (lastId == firstId) {
                    break;
                }
                triIds->InsertNextId(lastId);
            }
        }
    }
    return tri;
}

vtkUnstructuredGrid* buildLocalTangentMap(vtkPolyData* data) {
    vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::New();
    
    vtkNew<vtkCellArray> cells;
    for (size_t j = 0; j < data->GetNumberOfPoints(); j++) {
        vtkPolygon* tri = computeLocalTangentMap(data, j);
        cells->InsertNextCell(tri);
    }
    ugrid->SetCells(VTK_POLYGON, cells.GetPointer());
    return ugrid;
}


// Compute geodesic distance to the set of destination points
void computeGeodesicDistance(vtkPolyData* data, vtkDataArray* dest) {
    vGraph vg(data);
    vg.buildAdjacencyList();
    
    vGraph::vtkIdVector source;
    for (size_t u = 0; u < vg.adjList.size(); u++) {
        int ul = dest->GetTuple1(u);
        for (size_t j = 0; j < vg.adjList[u].size(); j++) {
            size_t v = vg.adjList[u][j];
            int vl = dest->GetTuple1(v);
            if (vl != ul) {
                source.push_back(u);
                break;
            }
        }
    }
    
    vg.computeGeodesicDistance(source);
    
    vtkDoubleArray* dist = vtkDoubleArray::New();
    dist->SetName("distance");
    dist->SetNumberOfComponents(1);
    dist->SetArray(&vg.dist[0], vg.dist.size(), 1);
    
    data->GetPointData()->AddArray(dist);
}


void processGeodesicOptions(pi::Options& opts) {
    opts.addOption("-computeGeodesicDistance", "Compute geodesic distance to the set of points", "-computeGeodesicDistance input-data -scalarName destinationScalar", SO_NONE);
    
    opts.addOption("-computeLocalExpMap", SO_NONE);
}

void processGeodesicCommands(pi::Options& opts, pi::StringVector& args) {
    if (opts.GetBool("-computeGeodesicDistance")) {
        string inputFile = args[0];
        string outputFile = args[1];
        
        string scalarName = opts.GetString("-scalarName", "BorderPoints");
        vtkPolyData* data = vio.readFile(inputFile);
        vtkDataArray* dest = data->GetPointData()->GetArray(scalarName.c_str());
        
        if (data == NULL || dest == NULL) {
            cout << "Input error!" << endl;
            return;
        }
        computeGeodesicDistance(data, dest);
        vio.writeFile(outputFile, data);
    } else if (opts.GetBool("-computeLocalExpMap")) {
        string inputFile = args[0];
        string outputFile = args[1] == "" ? "localmap.vtu" : args[1];

        vtkPolyData* data = vio.readFile(inputFile.c_str());
        if (data == NULL) {
            cout << "Input error: " << inputFile << endl;
            return;
        }
        vtkUnstructuredGrid* ugrid = buildLocalTangentMap(data);
        vio.writeFile(outputFile, ugrid);
    }
}