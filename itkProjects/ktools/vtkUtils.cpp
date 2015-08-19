//
//  vtkUtils.cpp
//  ktools
//
//  Created by Joowhi Lee on 8/17/15.
//
//

#include "vtkUtils.h"

#include <vtkPointData.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkCellLocator.h>
#include <vtkCellTreeLocator.h>
#include <vtkIntArray.h>
#include <vtkDoubleArray.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkNew.h>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include "vtkio.h"

using namespace std;

vtkIntArray* selectBoundaryCells(vtkDataSet* data, string scalarName) {
    const size_t nCells = data->GetNumberOfCells();
	vtkDataArray* array = data->GetPointData()->GetArray(scalarName.c_str());
	vector<int> scalarValue;
	vtkIntArray* intArray = vtkIntArray::New();
	intArray->SetNumberOfComponents(1);
	intArray->SetNumberOfTuples(nCells);
	
    for (size_t j = 0; j < nCells; j++) {
        vtkCell* cell = data->GetCell(j);
		scalarValue.clear();
        for (size_t k = 0; k < cell->GetNumberOfPoints(); k++) {
			vtkIdType cellPtId = cell->GetPointId(k);
			scalarValue.push_back(array->GetTuple1(cellPtId));
        }
		
		bool borderCell = false;
		if (scalarValue.size() > 0) {
			int value = scalarValue[0];
			for (size_t j = 1; j < scalarValue.size(); j++) {
				if (scalarValue[j] != value) {
					borderCell = true;
				}
			}
		}
		
		if (borderCell) {
			intArray->SetValue(j, 1);
		}
    }
	
	intArray->SetName("BorderCells");
	data->GetCellData()->AddArray(intArray);
	
	return intArray;
}


// mark the boundary points (2: interior boundary, 3: exterior boundary)
vtkDataArray* selectBoundaryPoints(vtkDataSet* data, std::string scalarName) {
    
    vtkDataArray* interiorMaker = data->GetPointData()->GetArray(scalarName.c_str());
    if (interiorMaker == NULL) {
        cout << "Can't find scalar values: " << scalarName << endl;
        return NULL;
    }
    
    const size_t npts = data->GetNumberOfPoints();
    
    vtkIntArray* boundaryMarker = vtkIntArray::New();
    boundaryMarker->SetNumberOfComponents(1);
    boundaryMarker->SetNumberOfTuples(npts);
	boundaryMarker->FillComponent(0, 0);
	

	unordered_set<vtkIdType> nbrs;
    vtkNew<vtkIdList> cellIds;

    for (size_t j = 0; j < npts; j++) {
        if (interiorMaker->GetTuple1(j) != 1) {
            continue;
        }
        
        // iterate over neighbor cells and find neighbor points
		nbrs.clear();
        cellIds->Reset();
        data->GetPointCells(j, cellIds.GetPointer());
        for (size_t k = 0; k < cellIds->GetNumberOfIds(); k++) {
            vtkCell* cell = data->GetCell(cellIds->GetId(k));
            for (size_t l = 0; l < cell->GetNumberOfEdges(); l++) {
                vtkCell* edge = cell->GetEdge(l);
                vtkIdType s = edge->GetPointId(0);
                vtkIdType e = edge->GetPointId(1);
                if (s == j) {
					nbrs.insert(e);
                } else if (e == j) {
					nbrs.insert(s);
                }
            }
        }
        
        // check neighbor points and find exterior points
        bool surfacePoint = false;
		unordered_set<vtkIdType>::const_iterator iter = nbrs.begin();
		for (; iter != nbrs.end(); iter++) {
			// if the neighbor is outside point
			vtkIdType nbrId = *iter;
            if (interiorMaker->GetTuple1(nbrId) == 0) {
                boundaryMarker->SetTuple1(nbrId, 3);
				// then j is a surface point
                surfacePoint = true;
            }
        }
        if (surfacePoint) {
            boundaryMarker->SetTuple1(j, 2);
        } else {
            boundaryMarker->SetTuple1(j, 1);
        }
    }
	
    boundaryMarker->SetName("BorderPoints");
    data->GetPointData()->AddArray(boundaryMarker);
    return boundaryMarker;
}


// using the scalar values, find edges in the grid that intersects the given surface
void selectSurfaceIntersectingEdges(vtkDataSet* grid, string scalarName, vtkPolyData* surf) {
    
}


void sampleSurfaceScalar(vtkDataSet* grid, string markerScalar, vtkPolyData* surf, string surfScalarName) {
    
    if (surf == NULL || grid == NULL) {
        cout << "sampleSurfaceScalar(): An input is NULL: " << grid << ", " << surf << endl;
        return;
    }
    vtkNew<vtkCellLocator> surfCellLoc;
    surfCellLoc->SetDataSet(surf);
    surfCellLoc->AutomaticOn();
    surfCellLoc->BuildLocator();

    vtkDataArray* markers = grid->GetPointData()->GetArray(markerScalar.c_str());
    vtkDataArray* surfScalars = surf->GetPointData()->GetArray(surfScalarName.c_str());
    unordered_set<int> cellScalars(10);
    
    const size_t nPts = markers->GetNumberOfTuples();
    
    vtkIntArray* outputScalar = vtkIntArray::New();
    outputScalar->SetNumberOfComponents(1);
    outputScalar->SetNumberOfValues(nPts);
    outputScalar->SetName("SampledSurfaceScalars");
    
    vtkIntArray* surfCheckedCells = vtkIntArray::New();
    surfCheckedCells->SetNumberOfComponents(1);
    surfCheckedCells->SetNumberOfValues(surf->GetNumberOfCells());
    surfCheckedCells->SetName("SampledSurfaceCells");
    
    for (size_t j = 0; j < nPts; j++) {
        // exterior border is 3
        int pointScalar = 0;
        
        // if the grid point j is an exterior boundary point
        if (markers->GetTuple1(j) == 3) {
            double pts[3];
            grid->GetPoint(j, pts);
            
            // find the closest cell from the grid point j (pts)
            double closestPoint[3] = { 0, };
            vtkIdType closestCellId = 0;
            int subId = -1;
            double dist2 = 0;
            
            // what if the closest point is a vertex?
            surfCellLoc->FindClosestPoint(pts, closestPoint, closestCellId, subId, dist2);
            
            // report this cell was checekd
            surfCheckedCells->SetValue(closestCellId, 1);
            
            
            // find the closest cell of the surface from the grid point j
            vtkCell* surfCell = surf->GetCell(closestCellId);
            if (surfCell == NULL) {
                cout << "Can't find closest cell for: " << j << endl;
            } else {
                // inspect the point scalar values of the closest cell
                const size_t nCellPts = surfCell->GetNumberOfPoints();
                int scalars[4] = {0,};
                for (size_t k = 0; k < nCellPts; k++) {
                    vtkIdType cellPointId = surfCell->GetPointId(k);
                    int surfScalar = (int) surfScalars->GetTuple1(cellPointId);
                    scalars[surfScalar]++;
                }
                for (size_t k = 0; k < 4; k++) {
                    if (scalars[k] >= 2) {
                        pointScalar = k ;
                        break;
                    }
                }
            }
            outputScalar->SetValue(j, pointScalar);
        } else if (markers->GetTuple1(j) == 2 || markers->GetTuple1(j) == 1) {
            outputScalar->SetValue(j, -1);
        }
    }
    
    grid->GetPointData()->AddArray(outputScalar);
    surf->GetCellData()->AddArray(surfCheckedCells);
    
}


void findNeighborPoints(vtkDataSet* data, vtkIdType ptId, unordered_set<vtkIdType>& nbrs, vtkIdList* cellIds) {
    cellIds->Reset();
    data->GetPointCells(ptId, cellIds);
    for (size_t k = 0; k < cellIds->GetNumberOfIds(); k++) {
        vtkCell* cell = data->GetCell(cellIds->GetId(k));
        for (size_t l = 0; l < cell->GetNumberOfEdges(); l++) {
            vtkCell* edge = cell->GetEdge(l);
            vtkIdType s = edge->GetPointId(0);
            vtkIdType e = edge->GetPointId(1);
            if (s == ptId) {
                nbrs.insert(e);
            } else if (e == ptId) {
                nbrs.insert(s);
            }
        }
    }
}

void buildAdjacencyList(vtkDataSet* surf, string roiLabel, vector<pair<vtkIdType,vector<vtkIdType> > >& graph) {
    
    // the roiArray extract the graph nodes within the ROI
    // Currently, I assume the roi is marked as -1 as I did in sampleSurfaceScalar
    vtkDataArray* roiArray = surf->GetPointData()->GetArray(roiLabel.c_str());
    graph.reserve(roiArray->GetNumberOfTuples());
    graph.clear();
    
    unordered_set<vtkIdType> nbrs;
    vtkNew<vtkIdList> cellIds;
    for (size_t j = 0; j < roiArray->GetNumberOfTuples(); j++) {
        if (roiArray->GetTuple1(j) == -1) {
            //  find the neighborhood of j
            findNeighborPoints(surf, j, nbrs, cellIds.GetPointer());
            graph.push_back(make_pair(j, vector<vtkIdType>(nbrs.begin(), nbrs.end())));
        }
        if (j % 1000 == 1) {
            cout << "Processed " << (j+1) << " points ..." << endl;
        }
    }
    
    for (size_t j = 0; j < graph.size(); j++) {
        cout << graph[j].first << " : ";
        for (size_t k = 0; k < graph[j].second.size(); k++) {
            cout << graph[j].second[k] << " ";
        }
        cout << endl;
    }
}


// select intersecting cells of 'surf' with 'grid' that has a scalar value 1 of 'scalarName'
vtkIntArray* selectIntersectingCells(vtkDataSet* surf, vtkDataSet* grid, std::string scalarName) {
	
	vtkNew<vtkCellTreeLocator> cellLoc;
	cellLoc->SetDataSet(surf);
	cellLoc->AutomaticOn();
	cellLoc->BuildLocator();
	
	vtkDataArray* scalars = grid->GetCellData()->GetArray(scalarName.c_str());
	
	for (size_t j = 0; j < scalars->GetNumberOfTuples(); j++) {
		vtkCell* cell = grid->GetCell(j);
        
        
	}
}




void processVTKUtilsOptions(pi::Options& opts) {
    opts.addOption("-markBorderCells", "Mark border cells of an input dataset. The border cells have 1 in BorderCells data", "-markBorderCells input-data output-data", SO_NONE);
    opts.addOption("-markBorderPoints", "Mark border points of an input dataset. The border points will be marked as 2 and its exterior neighbors will be marked as 3.", "-markBorderPoints input-data output-data", SO_NONE);
    opts.addOption("-sampleSurfaceScalars", "For each point marked as 3, sample the closest cell's majority scalar value.", "-sampleSurfaceScalars input-dataset input-surface output-dataset (-o output-surface)", SO_NONE);
    opts.addOption("-buildAdjacencyList", "Build an adjacency list for the input data and its point scalars", "-buildAdjacencyList input-data scalar-name", SO_NONE);
}


void processVTKUtils(pi::Options opts, pi::StringVector args) {
	vtkIO vio;
	string input1File, input2File, outputFile;
	if (opts.GetBool("-markBorderCells")) {
		input1File = args[0];
		outputFile = args[1];
		vtkDataSet* data1 = vio.readDataFile(input1File);
		string scalarName = opts.GetString("-scalarName", "SelectedPoints");
		selectBoundaryCells(data1, scalarName);
		vio.writeFile(outputFile, data1);
    } else if (opts.GetBool("-markBorderPoints")) {
		input1File = args[0];
		outputFile = args[1];
        vtkDataSet* data1 = vio.readDataFile(input1File);
        string scalarName = opts.GetString("-scalarName", "SelectedPoints");
        selectBoundaryPoints(data1, scalarName);
        vio.writeFile(outputFile, data1);
    } else if (opts.GetBool("-sampleSurfaceScalars")) {
        input1File = args[0];
        input2File = args[1];
        outputFile = args[2];
        vtkDataSet* data1 = vio.readDataFile(input1File);
        vtkPolyData* surf1 = vio.readFile(input2File);
        sampleSurfaceScalar(data1, opts.GetString("-scalarName", "BorderPoints").c_str(), surf1, "labels");
        vio.writeFile(outputFile, data1);
        if (opts.GetString("-o") != "") {
            vio.writeFile(opts.GetString("-o"), surf1);
        }
    } else if (opts.GetBool("-buildAdjacencyList")) {
        input1File = args[0];
        string scalarName = opts.GetString("-scalarName", "SampledSurfaceScalars");
        vtkDataSet* data1 = vio.readDataFile(input1File);
        
        vector<pair<vtkIdType, vector<vtkIdType> > > graph;
        buildAdjacencyList(data1, scalarName, graph);
    }
	
	
}
