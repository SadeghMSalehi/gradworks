#include "itkImageCommon.h"
#include "iostream"
#include "math.h"

using namespace std;

typedef itk::Image<double,3> ImageType;
typedef itk::ImageRegionIterator<ImageType> IteratorType;

int main(int argc, char* argv[]) {
	int ret;
	ImageType::Pointer srcImg = ReadImageT<ImageType>(argv[1], ret);
	IteratorType iter(srcImg, srcImg->GetRequestedRegion());

	for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
		double v = iter.Get();
		if (v < 0 || v > M_PI) {
			cout << v << endl;
		}
	}
}

