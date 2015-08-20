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



using namespace pi;
using namespace std;

static vtkIO vio;




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
    }
}