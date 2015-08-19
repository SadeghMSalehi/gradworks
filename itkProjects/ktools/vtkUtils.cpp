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
#include <vtkNew.h>
#include <vector>

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
	}
	
	
}
