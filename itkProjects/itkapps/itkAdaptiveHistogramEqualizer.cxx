#include "itkImageCommon.h"
#include "iostream"
#include "itkAdaptiveHistogramEqualizationImageFilter.h"

using namespace std;

typedef itk::Image<double,3> ImageType;
typedef itk::AdaptiveHistogramEqualizationImageFilter<ImageType> EqualiziationFilterType;

int main(int argc, char* argv[]) {
	if (argc < 6) {
		cout << "usage: " << argv[0] << " input-image output-image alpha beta window" << endl;
		return 0;
	}

	int ret = 0;
	ImageType::Pointer imgPtr = ReadImageT<ImageType>(argv[1], ret);

	EqualiziationFilterType::Pointer filter = EqualiziationFilterType::New();
	filter->SetInput(imgPtr);
	filter->SetAlpha(atof(argv[3]));
	filter->SetBeta(atof(argv[4]));

	ImageType::SizeType r;
	r[0] = atoi(argv[5]);
	r[1] = atoi(argv[5]);
	r[2] = atoi(argv[5]);
	filter->SetRadius(r);
	filter->Update();

	WriteImageT<ImageType>(argv[2], filter->GetOutput());
}
