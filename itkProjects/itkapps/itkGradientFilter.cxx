#include "itkImageCommon.h"
#include "itkGradientImageFilter.h"
#include "iostream"

using namespace std;

typedef itk::Image<double,3> ImageType;
typedef itk::GradientImageFilter<ImageType,double,double> GradientFilterType;

int main(int argc, char* argv[]) {
	if (argc < 3) {
		cout << "usage: " << argv[0] << " input-image output-image " << endl;
		return 0;
	}

	int ret;
	ImageType::Pointer srcImg = ReadImageT<ImageType>(argv[1], ret);
	GradientFilterType::Pointer gf = GradientFilterType::New();
	gf->SetInput(srcImg);
	gf->Update();
	WriteImageT<GradientFilterType::OutputImageType>(argv[2], gf->GetOutput());
}
