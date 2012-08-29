#include "itkImageCommon.h"
#include "iostream"

typedef itk::Image<int, 3> IntImage;

using namespace std;

int main(int argc, char* argv[]) {
	int ret = 0;
	IntImage::Pointer img = ReadImageT<IntImage>(argv[1], ret);
	int* buff = img->GetBufferPointer();
	cout << buff[0] << endl;
}
