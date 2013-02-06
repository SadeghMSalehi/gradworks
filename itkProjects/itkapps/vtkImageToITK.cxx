#include "itkImageCommon.h"
#include "vtkXMLImageDataReader.h"
#include "vtkImageData.h"
#include "itkImage.h"

typedef itk::Image<double,3> ImageType;

int main(int argc, char* argv[]) {
	vtkXMLImageDataReader* r = vtkXMLImageDataReader::New();
	r->SetFileName(argv[1]);
	r->Update();

	vtkImageData* vti = r->GetOutput();
	int* dims = vti->GetDimensions();

	ImageType::Pointer newImg = NewImageT<ImageType>(dims[0], dims[1], dims[2]);

	for (int z = 0; z < dims[2]; z++) {
		for (int y = 0; y < dims[1]; y++) {
			for (int x = 0; x < dims[0]; x++) {
				double v = vti->GetScalarComponentAsDouble(x,y,z,0);
				ImageType::IndexType idx;
				idx[0] = x;
				idx[1] = y;
				idx[2] = z;
				newImg->SetPixel(idx, v);
			}
		}
	}

	WriteImageT<ImageType>(argv[2], newImg);
}
