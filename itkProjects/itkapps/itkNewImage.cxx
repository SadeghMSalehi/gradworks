#include "itkImageCommon.h"

typedef itk::Image<double,3> ImageType;

using namespace std;

int main(int argc, char* argv[]) {
	if (argc < 4) {
		cout << "Usage: " << argv[0] << " x-size y-size z-size" << endl;
		return 0;
	}

	ImageType::Pointer img = NewImageT<ImageType>(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
	WriteImageT<ImageType>("tempImage.nrrd", img);	
}
