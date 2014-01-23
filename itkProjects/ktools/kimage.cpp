//
//  kimage.cpp
//  ktools
//
//  Created by Joohwi Lee on 1/22/14.
//
//

#include "kimage.h"

void SetArrayTuple(vtkDataArray* a, int i, VectorType v) {
    for (int j = 0; j < a->GetNumberOfComponents(); j++) {
        a->SetComponent(i, j, v[j]);
    }
}

void SetArrayTuple(vtkDataArray* a, int i, double v) {
    a->SetTuple1(i, v);
}

