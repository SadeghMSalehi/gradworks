#include <iostream>
#include <ostream>
#include <itkResampleImageFilter.h>
#include <itkAffineTransform.h>
#include <itkMath.h>
#include "itkImageCommon.h"

using namespace std;
using namespace itk;
using namespace itkcmds;

typedef itk::Image<double,3> ImageType;
typedef itk::AffineTransform<double,3> TransformType;
typedef itk::ResampleImageFilter<ImageType, ImageType> ResamplerType;

int main(int argc, char* argv[]) {
	itkImageIO<ImageType> imageIO;
	ImageType::Pointer srcImg = imageIO.ReadImageT(argv[1]);

	TransformType::Pointer txfm = TransformType::New();
	ResamplerType::LinearInterpolatorPointerType linIntp = ResamplerType::LinearInterpolatorType::New();
	ResamplerType::Pointer resampler = ResamplerType::New();

	ImageType::SizeType srcSize = srcImg->GetLargestPossibleRegion().GetSize();
	ImageType::IndexType centerIdx;
	centerIdx[0] = srcSize[0] / 2;
	centerIdx[1] = srcSize[1] / 2;
	centerIdx[2] = srcSize[2] / 2;
	ImageType::PointType centerPoint;
	srcImg->TransformIndexToPhysicalPoint(centerIdx, centerPoint);

	TransformType::OutputVectorType txyz;
	txyz[0] = 1;
	txyz[1] = 0;
	txyz[2] = 0;
	//txfm->Translate(txyz);

	txfm->Rotate(0, 1, 10. * itk::Math::pi / 180.);

	resampler->SetInput(srcImg);
	resampler->SetTransform(static_cast<ResamplerType::TransformType::Pointer>(txfm));
	resampler->SetInterpolator(linIntp);
	resampler->SetReferenceImage(srcImg);
	resampler->UseReferenceImageOn();
	resampler->Update();

	ImageType::Pointer outImg = resampler->GetOutput();
	cout << "Image Output Size: " << outImg->GetLargestPossibleRegion().GetSize() << endl;
	imageIO.WriteImageT(argv[2], outImg);
	
}
