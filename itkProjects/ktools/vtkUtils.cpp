//
//  vtkUtils.cpp
//  ktools
//
//  Created by Joowhi Lee on 8/17/15.
//
//

#include "vtkUtils.h"

void selectBoundaryCells(vtkDataSet* data) {
    const size_t nCells = data->GetNumberOfCells();
    for (size_t j = 0; j < nCells; j++) {
        vtkCell* cell = data->GetCell(j);
        for (size_t k = 0; k < cell->GetNumberOfPoints(); k++) {
            
        }
    }
}