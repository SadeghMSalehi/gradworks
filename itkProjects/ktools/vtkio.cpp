//
//  vtkio.cpp
//  ktools
//
//  Created by Joohwi Lee on 12/5/13.
//
//

#include "vtkio.h"
#include <vtkDataSet.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>

static bool endswith(std::string file, std::string ext) {
    int epos = file.length() - ext.length();
    if (epos < 0) {
        return false;
    }
    return file.rfind(ext) == epos;
}


/// @brief Read a vtk/vtp file. The file's type is automatically determined by its extension.
vtkPolyData* vtkIO::readFile(std::string file) {
    if (endswith(file, ".vtp")) {
        vtkXMLPolyDataReader* r = vtkXMLPolyDataReader::New();
        r->SetFileName(file.c_str());
        r->Update();
        return r->GetOutput();
    } else if (endswith(file, ".vtk")) {
        vtkPolyDataReader* r = vtkPolyDataReader::New();
        r->SetFileName(file.c_str());
        r->Update();
        return r->GetOutput();
    }
}

/// @brief Write a vtk/vtp file. The file's type is automatically determined by its extension.
void vtkIO::writeFile(std::string file, vtkDataSet *mesh) {
    if (endswith(file, ".vtp")) {
        vtkXMLPolyDataWriter* w = vtkXMLPolyDataWriter::New();
        w->SetInput(mesh);
        w->SetFileName(file.c_str());
        w->Write();
        w->Delete();
    } else if (endswith(file, ".vtk")) {
        vtkPolyDataWriter* w = vtkPolyDataWriter::New();
        w->SetInput(mesh);
        w->SetFileName(file.c_str());
        w->Write();
        w->Delete();
    } else if (endswith(file, ".vtu")) {
        vtkXMLUnstructuredGridWriter* w = vtkXMLUnstructuredGridWriter::New();
        w->SetInput(mesh);
        w->SetFileName(file.c_str());
        w->Write();
        w->Delete();
    }
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
