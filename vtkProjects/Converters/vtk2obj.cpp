#include "vtkPolyDataReader.h"
#include "vtkMNIObjectWriter.h"

int main(int argc, char* argv[]) {
	if (argc < 3) {
		printf("usage: %s input.vtk output.obj\n", argv[0]);
		return 0;
	}
	vtkPolyDataReader* r = vtkPolyDataReader::New();
	r->SetFileName(argv[1]);
	r->Update();
	vtkMNIObjectWriter* w = vtkMNIObjectWriter::New();
	w->SetFileName(argv[2]);
	w->SetInputConnection(r->GetOutputPort());
	w->Write();
}
