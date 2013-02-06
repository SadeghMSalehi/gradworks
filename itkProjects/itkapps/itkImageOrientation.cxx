#include "itkImageCommon.h"
#include "itkOrientedImage.h"
#include "iostream"
#include "math.h"

using namespace std;

typedef itk::Image<unsigned short,3> ImageType;
typedef itk::OrientedImage<unsigned short, 3> OrientedImageType;

int main(int argc, char* argv[]) {
	int ret;
	ImageType::Pointer srcImg = ReadImageT<ImageType>(argv[1], ret);
	cout << srcImg->GetDirection() << endl;
	OrientedImageType::Pointer srcImg2 = ReadImageT<OrientedImageType>(argv[1], ret);
	cout << srcImg2->GetDirection() << endl;
}

