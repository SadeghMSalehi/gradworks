#include "itkImageCommon.h"
#include "itkGradientMagnitudeImageFilter.h"
#include "iostream"

using namespace std;

typedef itk::Image<unsigned short, 3> ImageType;
typedef itk::Image<double, 3> OutImageType;

typedef itk::GradientMagnitudeImageFilter<ImageType,OutImageType> GradientMagnitudeFilter;

int main(int argc, char* argv[]) {
	if (argc < 3) {
		cout << "usage: " << argv[0] << " input-image output-image" << endl;
		return 0;
	}

	int ret;
	ImageType::Pointer inputImg;
	OutImageType::Pointer outImg;

	inputImg = ReadImageT<ImageType>(argv[1], ret);

	GradientMagnitudeFilter::Pointer gmf = GradientMagnitudeFilter::New();
	gmf->SetInput(inputImg);
	gmf->Update();
	outImg = gmf->GetOutput();

	WriteImageT<OutImageType>(argv[2], outImg);
}
