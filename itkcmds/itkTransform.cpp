#include <iostream>
#include <itkLinearInterpolateImageFunction.h>
#include <itkResampleImageFilter.h>
#include <itkAffineTransform.h>
#include "itkImageCommon.h"

using namespace std;
using namespace itk;
using namespace itkcmds;

typedef itk::Image<double,3> ImageType;
typedef itk::AffineTransform<double, 3> TransformType;
typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;
typedef itk::ResampleImageFilter<ImageType, ImageType> ResamplerType;

int main(int argc, char* argv[]) {
	itkImageIO<ImageType> imageIO;
	ImageType::Pointer srcImg = imageIO.ReadImageT(argv[1]);

	cout << "Image information" << endl;
	cout << "Origin: " << srcImg->GetOrigin() << endl;
	cout << "Spacing: " << srcImg->GetSpacing() << endl;
	cout << "Size: " << srcImg->GetLargestPossibleRegion().GetSize() << endl;

	
	TransformType::Pointer txfm = TransformType::New();
	InterpolatorType::Pointer intp = InterpolatorType::New();
	ResamplerType::Pointer resampler = ResamplerType::New();

	TransformType::InputPointType center;
	center[0] = 0;
	center[1] = 0;
	center[2] = 0;
	txfm->SetCenter(center);

	TransformType::OutputVectorType txyz;
	txyz[0] = 10;
	txyz[1] = 0;
	txyz[2] = 0;
	txfm->Translate(txyz);

	TransformType::OutputVectorType scales;

	cout << "Translation: " << txyz << endl;
	cout << "Scales: " << scales << endl;
	cout << "Transform: " << txfm << endl;

	resampler->SetInput(srcImg);
	resampler->SetTransform(txfm);
	resampler->SetInterpolator(intp);
	resampler->Update();

	ImageType::Pointer outImg = resampler->GetOutput();
	//imageIO.WriteImageT(argv[2], outImg);
	
}