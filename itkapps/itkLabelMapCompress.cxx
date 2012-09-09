#include "itkImageCommon.h"
#include "iostream"

using namespace std;

int main(int argc, char* argv[]) {
  typedef itk::Image<unsigned int,3> ImageType;
	if (argc < 3) {
		cout << "usage: " << argv[1] << " input-image output-image" << endl;
		return 0;
	}
	int ret = 0;
	WriteImageT<ImageType>(argv[2], ReadImageT<ImageType>(argv[1], ret));
	return 0;
}
