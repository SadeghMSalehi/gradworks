#include <iostream>
#include "itkImageCommon.h"

typedef itk::Image<int, 3> IntImage;
typedef itk::Image<float, 3> FloatImage;

using namespace std;
using namespace itkcmds;
using namespace itk;

int main(int argc, char* argv[]) {
	cout << "itkinfo version 0.1" << endl;
	if (argc < 1) {
		cout << "usage: itkinfo filename" << endl;
		return 0;
	}

	itkImageIO<FloatImage> imageFloatIO;
	itkImageIO<IntImage> imageIntIO;

	FloatImage::Pointer img = imageFloatIO.ReadImageT(argv[1]);
	if (img.IsNull()) {
		cout << "Failed in loading [" << argv[1] << "]" << endl;
	} else {
		cout << "Printing Image Information" << endl;
		img->Print(cout);
		cout << endl;
	}


	return 0;
}
