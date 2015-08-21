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


void computeLocalTangentMap(vtkPolyData* g, vtkUnstructuredGrid* outputTangentMaps) {
    const size_t nPts = g->GetNumberOfPoints();

    // output data
    vtkIdTypeArray* ringIds = vtkIdTypeArray::New();
    ringIds->SetName("RingIDs");
    ringIds->SetNumberOfComponents(1);
    
    vtkCellArray* tangentMapArray = vtkCellArray::New();
    tangentMapArray->Allocate(nPts);
    
    vtkPoints* tangentPoints = vtkPoints::New();
    tangentPoints->Allocate(nPts*6);
    
    
    // temporary data types
    vtkNew<vtkIdList> cellIds, ptIds;
    vector<deque<vtkIdType> > neighbors;
    for (size_t p = 0; p < nPts; p++) {
        cellIds->Reset();
        ptIds->Reset();
        
        
        // inspect the neighbor cells of point p
        g->GetPointCells(p, cellIds.GetPointer());
        neighbors.resize(cellIds->GetNumberOfIds());
        
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
        
        // construct a tangent plane with neighbor points
        vtkPolygon* tangentPlane = vtkPolygon::New();
        
        vtkIdType firstId = neighbors[0][1];
        vtkIdType lastId = -1;
        for (size_t j = 1; j < neighbors[0].size(); j++) {
            lastId = neighbors[0][j];
            ringIds->InsertNextValue(lastId);
            vtkIdType pId = tangentPoints->InsertNextPoint(g->GetPoint(lastId));
            tangentPlane->GetPointIds()->InsertNextId(pId);
        }
        
        for (size_t j = 0; j < neighbors.size() && lastId != firstId; j++) {
            for (size_t k = 0; k < neighbors.size(); k++) {
                if (neighbors[k][1] == lastId) {
                    lastId = neighbors[k][2];
                    if (lastId == firstId) {
                        break;
                    }
                    ringIds->InsertNextValue(lastId);
                    vtkIdType pId = tangentPoints->InsertNextPoint(g->GetPoint(lastId));
                    tangentPlane->GetPointIds()->InsertNextId(pId);
                }
            }
        }
        
        if (p % 100 == 0) {
            cout << "Points processed: " << p << " ..." << endl;
        }
    }

    outputTangentMaps->SetCells(VTK_POLYGON, tangentMapArray);
    outputTangentMaps->GetCellData()->AddArray(ringIds);
    outputTangentMaps->SetPoints(tangentPoints);
    
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
    
    opts.addOption("-computeLocalTangentMap", SO_NONE);
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
    } else if (opts.GetBool("-computeLocalTangentMap")) {
        string inputFile = args[0];
        string outputFile = args[1] == "" ? "localmap.vtu" : args[1];

        vtkPolyData* data = vio.readFile(inputFile.c_str());
        if (data == NULL) {
            cout << "Input error: " << inputFile << endl;
            return;
        }
        vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::New();
        computeLocalTangentMap(data, ugrid);
        vio.writeFile(outputFile, ugrid);
    }
}