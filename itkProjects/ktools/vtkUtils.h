//
//  vtkUtils.h
//  ktools
//
//  Created by Joowhi Lee on 8/17/15.
//
//

#ifndef __ktools__vtkUtils__
#define __ktools__vtkUtils__

#include <stdio.h>
#include <string>


#include "piOptions.h"


class vtkIntArray;
class vtkPolyData;
class vtkDataArray;
class vtkDataSet;

vtkIntArray* selectBoundaryCells(vtkDataSet* data, std::string scalarName);
vtkDataArray* selectBoundaryPoints(vtkDataSet* data, std::string scalarName);
vtkIntArray* selectIntersectingCells(vtkDataSet* surf, vtkDataSet* grid, std::string scalarName);


void processVTKUtilsOptions(pi::Options& opts);

void processVTKUtils(pi::Options opts, pi::StringVector args);

#endif /* defined(__ktools__vtkUtils__) */
