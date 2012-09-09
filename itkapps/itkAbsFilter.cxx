#include "itkImageCommon.h"
#include "itkAbsImageFilter.h"

using namespace std;

typedef itk::Image<double,3> ImageType;
typedef itk::AbsImageFilter<ImageType,ImageType> AbsFilterType;

int main(int argc, char* argv[]) {
	if (argc < 3) {
		cout << "usage: " << argv[0] << " input-image output-image" << endl;
		return 0;
	}
	
	int ret = 0;
	ImageType::Pointer srcImg = ReadImageT<ImageType>(argv[1], ret);
	AbsFilterType::Pointer filt = AbsFilterType::New();
	filt->SetInput(srcImg);
	filt->Update();
	WriteImageT<ImageType>(argv[2], filt->GetOutput());
}
