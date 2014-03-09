#include "itkImageCommon.h"
#include "iostream"

using namespace std;

typedef itk::Image<double,3> ImageType;

int main(int argc, char* argv[]) {
	if (argc < 3) {
		cout << argv[0] << " input-image output-image" << endl;
		return 0;
	}

	int ret = 0;
	ImageType::Pointer img = ReadImageT<ImageType>(argv[1], ret);	
	ImageType::SizeType imgSz = img->GetBufferedRegion().GetSize();
	for (unsigned int i = 0; i < imgSz[2]; i++) {
		for (unsigned int j = 0; j < imgSz[1]; j++) {
			for (unsigned int k = 0; k < imgSz[0]; k++) {
				ImageType::IndexType idx;
				idx[0] = k;
				idx[1] = j;
				idx[2] = i;
				double x = img->GetPixel(idx);
				img->SetPixel(idx, 1-x);
			}
		}
	}
	WriteImageT<ImageType>(argv[2], img);
	return 0;
}
