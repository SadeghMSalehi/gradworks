#include "itkImageCommon.h"
#include "EdgePreprocessingImageFilter.h"

#include "iostream"

using namespace std;

typedef itk::Image<unsigned short, 3> ImageType;
typedef itk::Image<double, 3> OutImageType;
typedef EdgePreprocessingImageFilter<ImageType, OutImageType> EdgePreprocessingFilter;

int main(int argc, char* argv[]) {
	if (argc < 3) {
		cout << "usage: " << argv[0] << " input-image output-image" << endl;
		return 0;
	}

	int ret = 0;
	ImageType::Pointer inImage = ReadImageT<ImageType>(argv[1], ret);

	EdgePreprocessingSettings settings;
	settings.SetGaussianBlurScale(1.0f);
	settings.SetRemappingSteepness(0.05f);
	settings.SetRemappingExponent(3.0f);

	EdgePreprocessingFilter::Pointer filter = EdgePreprocessingFilter::New();
	filter->SetEdgePreprocessingSettings(settings);
	filter->SetInput(inImage);
	filter->Update();
	OutImageType::Pointer outImage = filter->GetOutput();

	WriteImageT<OutImageType>(argv[2], outImage);
}
