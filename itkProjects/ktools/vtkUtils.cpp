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
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkNew.h>
#include <vector>
#include <unordered_set>
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


// mark surface boundary and find intersecting cells
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
//            boundaryMarker->SetTuple1(j, 1);
        }
    }
	
    boundaryMarker->SetName("BorderPoints");
    data->GetPointData()->AddArray(boundaryMarker);
    return boundaryMarker;
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
}


void processVTKUtils(pi::Options opts, pi::StringVector args) {
	vtkIO vio;
	string input1File, input2File, outputFile;
	if (opts.GetBool("-markBorderCells")) {
		input1File = args[0];
		outputFile = args[1];
		vtkDataSet* data1 = vio.readDataFile(input1File);
		string scalarName = opts.GetString("-scalarName", "InteriorPoints");
		selectBoundaryCells(data1, scalarName);
		vio.writeFile(outputFile, data1);
    } else if (opts.GetBool("-markBorderPoints")) {
		input1File = args[0];
		outputFile = args[1];
        vtkDataSet* data1 = vio.readDataFile(input1File);
        string scalarName = opts.GetString("-scalarName", "InteriorPoints");
        selectBoundaryPoints(data1, scalarName);
        vio.writeFile(outputFile, data1);
    }
	
	
}
