#include "vtkMNIObjectReader.h"
#include "vtkPolyDataWriter.h"

int main(int argc, char* argv[]) {
	if (argc < 3) {
		printf("usage: %s input.obj output.vtk\n", argv[0]);
		return 0;
	}
	vtkMNIObjectReader* r = vtkMNIObjectReader::New();
	r->SetFileName(argv[1]);
	r->Update();
	vtkPolyDataWriter* w = vtkPolyDataWriter::New();
	w->SetFileName(argv[2]);
	w->SetInput(r->GetOutputPort());
	w->Write();
}
