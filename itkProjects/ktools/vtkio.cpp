//
//  vtkio.cpp
//  ktools
//
//  Created by Joohwi Lee on 12/5/13.
//
//

#include "vtkio.h"
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkXMLPolyDataWriter.h>

vtkPolyData* vtkIO::readFile(std::string file) {
    vtkPolyDataReader* r = vtkPolyDataReader::New();
    r->SetFileName(file.c_str());
    r->Update();
    return r->GetOutput();
}

void vtkIO::writeFile(std::string file, vtkPolyData *mesh) {
    vtkPolyDataWriter* w = vtkPolyDataWriter::New();
    w->SetInput(mesh);
    w->SetFileName(file.c_str());
    w->Write();
    w->Delete();
}


void vtkIO::writeXMLFile(std::string file, vtkPolyData *mesh) {
    vtkXMLPolyDataWriter* w = vtkXMLPolyDataWriter::New();
    w->SetInput(mesh);
    w->SetFileName(file.c_str());
    w->SetDataModeToAppended();
    w->EncodeAppendedDataOff();
    w->SetCompressorTypeToZLib();
    w->SetDataModeToBinary();
    w->Write();
    w->Delete();
}