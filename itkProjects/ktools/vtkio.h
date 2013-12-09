//
//  vtkio.h
//  ktools
//
//  Created by Joohwi Lee on 12/5/13.
//
//

#ifndef __ktools__vtkio__
#define __ktools__vtkio__

#include <iostream>

class vtkPolyData;

class vtkIO {
public:
    vtkPolyData* readFile(std::string file);
    void writeFile(std::string file, vtkPolyData* mesh);
};

#endif /* defined(__ktools__vtkio__) */