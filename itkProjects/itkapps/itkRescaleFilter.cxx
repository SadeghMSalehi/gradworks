#include "itkImageCommon.h"
#include "iostream"
#include "itkRescaleIntensityImageFilter.h"

using namespace std;

typedef itk::Image<double,3> ImageType;
typedef itk::RescaleIntensityImageFilter<ImageType,ImageType> RescaleFilterType;

int main(int argc, char* argv[]) {
	if (argc < 5) {
		cout << "usage: " << argv[0] << "input-image output-image min max" << endl;
		return 0;
	}

	int ret = 0;
	ImageType::Pointer imgPtr = ReadImageT<ImageType>(argv[1], ret);

	RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
	rescaleFilter->SetInput(imgPtr);
	rescaleFilter->SetOutputMinimum(atof(argv[3]));
	rescaleFilter->SetOutputMaximum(atof(argv[4]));
	rescaleFilter->Update();

	WriteImageT<ImageType>(argv[2], rescaleFilter->GetOutput());
}
