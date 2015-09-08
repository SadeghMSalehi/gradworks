//
//  vGraph.cpp
//  ktools
//
//  Created by Joowhi Lee on 8/20/15.
//
//

#include "vGraph.h"

#include <deque>

#include <vtkMath.h>

using namespace std;


static void findNeighborPoints(vtkDataSet* data, vtkIdType ptId, unordered_set<vtkIdType>& nbrs, vtkIdList* cellIds, vtkIdTypeArray* neighbors = NULL) {
    nbrs.clear();
    cellIds->Reset();
    data->GetPointCells(ptId, cellIds);
    for (size_t k = 0; k < cellIds->GetNumberOfIds(); k++) {
        vtkCell* cell = data->GetCell(cellIds->GetId(k));
        for (size_t l = 0; l < cell->GetNumberOfEdges(); l++) {
            vtkCell* edge = cell->GetEdge(l);
            vtkIdType s = edge->GetPointId(0);
            vtkIdType e = edge->GetPointId(1);
            vtkIdType n = -1;
            if (s == ptId) {
                n = e;
            } else if (e == ptId) {
                n = s;
            }
            if (n > -1) {
                if (nbrs.find(n) == nbrs.end()) {
                    nbrs.insert(n);
                    if (neighbors != NULL) {
                        neighbors->SetComponent(ptId, nbrs.size() - 1, n);
                    }
                }
            }
        }
    }
}


// build the adjacency list of the dataset
void vGraph::buildAdjacencyList() {
    const size_t nPoints = data->GetNumberOfPoints();
    adjList.resize(nPoints);
    
    unordered_set<vtkIdType> nbrs;
    vtkNew<vtkIdList> cellIds;
    
    for (vtkIdType j = 0; j < nPoints; j++) {
        findNeighborPoints(data, j, nbrs, cellIds.GetPointer());
        adjList[j].assign(nbrs.begin(), nbrs.end());
    }
}


// perform BSF-based geodesic distance computation
// make sure whether the adjacency list is already prepared
void vGraph::computeGeodesicDistance(vtkIdVector &source) {
    vector<Color> color(adjList.size());
    deque<vtkIdType> q;

    // initialization
    dist.resize(data->GetNumberOfPoints());
    pred.resize(data->GetNumberOfPoints());
    
    fill(dist.begin(), dist.end(), DBL_MAX);
    fill(pred.begin(), pred.end(), -1);
    
    // source point
    for (size_t j = 0; j < source.size(); j++) {
        // nodes in the queue is Gray
        vtkIdType s = source[j];
        color[s] = White;
        pred[s] = s;
        dist[s] = 0;
        q.push_back(s);
        cout << j << "," << q.size() << endl;
    }

    cout << q.size() << endl;
    while (!q.empty()) {
        vtkIdType u = q.front();
        q.pop_front();
        
        // the white nodes already have the shortest path
        color[u] = White;
        double udist = dist[u];
        
        double uPt[3];
        data->GetPoint(u, uPt);
        
        // perform BSF
        vtkIdVector& nbrs = adjList[u];
        vtkIdVector::const_iterator iter = nbrs.begin();
        for (; iter != nbrs.end(); iter++) {
            const vtkIdType v = *iter;
            if (color[v] == White) {
                continue;
            }
            
            double vdist = dist[v];
            double vPt[3];
            data->GetPoint(v, vPt);
            
            // compute the distance
			const double uvdist = sqrt(vtkMath::Distance2BetweenPoints(uPt, vPt));
            
            // relaxation
            if (udist + uvdist < vdist) {
                pred[v] = u;
                dist[v] = udist + uvdist;
            }
            
            if (color[v] == Black) {
                q.push_back(v);
                color[v] = Gray;
            }
        }
    }
    
}