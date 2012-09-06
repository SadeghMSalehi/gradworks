/*
 * itkImageReg.cpp
 *
 *  Created on: Sep 5, 2012
 *      Author: joohwi
 */

#include "itkImageCommon.h"
#include "itkExceptionObject.h"
#include "itkAffineTransform.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkImageRegistrationMethod.h"
#include "itkResampleImageFilter.h"
#include "itkContinuousIndex.h"
#include "itkMath.h"
#include "itkTimer.h"
#include "iostream"

#include "MatrixCode.h"

using namespace std;

typedef itk::ExceptionObject itkException;
typedef itk::Image<double, 3> ImageType;
typedef itk::AffineTransform<double, 3> TransformType;
typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
typedef itk::MeanSquaresImageToImageMetric<ImageType, ImageType> MetricType;
typedef itk::LinearInterpolateImageFunction<ImageType, double> InterpolatorType;
typedef itk::ImageRegistrationMethod<ImageType, ImageType> RegistrationType;
typedef itk::ResampleImageFilter<ImageType, ImageType> ResamplerType;

#define angle2rad(a) (a*itk::Math::pi/180.0)

void createMat4FromMat(const ImageType::DirectionType &mat,
		MathCode::Mat4& matOut) {
	for (int r = 0; r < 3; r++) {
		for (int c = 0; c < 3; c++) {
			matOut[r][c] = mat[r][c];
		}
	}
}

void createMat4FromVec(const ImageType::SpacingType &spacing,
		MathCode::Mat4& matOut, bool inverse = false) {
	for (int i = 0; i < 3; i++) {
		matOut[i][i] = spacing[i];
		if (inverse) {
			matOut[i][i] = 1 / matOut[i][i];
		}
	}
}

void createIndex2IndexTransformMatrix(ImageType::Pointer src,
		ImageType::Pointer ref, TransformType::MatrixType &transformMatrix,
		MathCode::Mat4& matOut) {
	MathCode::Mat4 srcDirection, srcSpacing, srcImageToWorld;
	createMat4FromMat(src->GetDirection(), srcDirection);
	createMat4FromVec(src->GetSpacing(), srcSpacing);
	MathCode::mult(srcDirection, srcSpacing, srcImageToWorld);
	for (int i = 0; i < 3; i++) {
		srcImageToWorld[i][3] = src->GetOrigin()[i];
	}

	MathCode::Mat4 srcTransform, srcToRefWorld;
	createMat4FromMat(static_cast<ImageType::DirectionType>(transformMatrix),
			srcTransform);
	mult(srcTransform, srcImageToWorld, srcToRefWorld);

	MathCode::Mat4 refDirection, refSpacing, refSpacingDirection,
			refWorldToImage;
	createMat4FromVec(ref->GetSpacing(), refSpacing, true);
	createMat4FromMat(ref->GetInverseDirection(), refDirection);

	mult(refSpacing, refDirection, refSpacingDirection);
	for (int i = 0; i < 3; i++) {
		srcToRefWorld[i][3] -= ref->GetOrigin()[i];
	}
	cout << "Src To Ref World: " << srcToRefWorld << endl;
	cout << "Ref Spacing: " << refSpacing << endl;
	cout << "Ref SpacingDirection: " << refSpacingDirection << endl;
	mult(refSpacingDirection, srcToRefWorld, refWorldToImage);
	matOut.copyFrom(refWorldToImage);
}

ImageType::Pointer resampleImageRaw(ImageType::Pointer src,
		ImageType::Pointer ref, TransformType::Pointer transform) {
	MathCode::Mat4 srcIndexToRefIndex, refIndexToSrcIndex;
	TransformType::MatrixType transformMatrix = transform->GetMatrix();
	createIndex2IndexTransformMatrix(src, ref, transformMatrix,
			srcIndexToRefIndex);

	ImageType::RegionType srcRegion = src->GetBufferedRegion();
	double* srcBuff = src->GetBufferPointer();
	ImageType::RegionType refRegion = ref->GetBufferedRegion();
	ImageType::SizeType refSize = refRegion.GetSize();

	itkcmds::itkImageIO<ImageType> imageIO;
	ImageType::Pointer dst = imageIO.NewImageT(src);
	dst->FillBuffer(0);
	double* dstBuff = dst->GetBufferPointer();
	srcIndexToRefIndex.inverse(refIndexToSrcIndex);

	InterpolatorType::Pointer interp = InterpolatorType::New();
	interp->SetInputImage(src);

	for (int z = 0; z < refSize[2]; z++) {
		for (int y = 0; y < refSize[1]; y++) {
			for (int x = 0; x < refSize[0]; x++) {
				MathCode::Vec4 dstPoint(x, y, z, 1);
				MathCode::Vec4 srcPoint;
				mult(refIndexToSrcIndex, dstPoint, srcPoint);
				itk::ContinuousIndex<double,3> srcIndex;
				srcIndex[0] = srcPoint[0];
				srcIndex[1] = srcPoint[1];
				srcIndex[2] = srcPoint[2];
				if (!srcRegion.IsInside(srcIndex)) {
					continue;
				}

			}
		}
	}

	return src;
}

ImageType::Pointer resampleImageITK(ImageType::Pointer src,
		TransformType::Pointer transform) {
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	ResamplerType::Pointer resampler = ResamplerType::New();
	resampler->SetInput(src);
	resampler->SetTransform(
			static_cast<ResamplerType::TransformType::Pointer>(transform));
	resampler->SetInterpolator(interpolator);
	resampler->SetReferenceImage(src);
	resampler->UseReferenceImageOn();
	resampler->Update();

	return resampler->GetOutput();
}

void computeMSE(ImageType::Pointer src, ImageType::Pointer dst) {
	double* srcBuff = src->GetBufferPointer();
	double* dstBuff = dst->GetBufferPointer();

	ImageType::RegionType srcRegion = src->GetBufferedRegion();
	int nPixels = srcRegion.GetNumberOfPixels();

	double mse = 0;
	int nCount = 0;
	for (int i = 0; i < nPixels; i++) {
		double diff = srcBuff[i] - dstBuff[i];
		mse += diff * diff;
		nCount++;
	}
	mse /= nPixels;
	cout << "Mean Squared Error: " << mse << "; # pixels: " << nCount << endl;
}

int main(int argc, char* argv[]) {
	itkcmds::itkImageIO<ImageType> imageIO;
	ImageType::Pointer srcImg = imageIO.ReadImageT(argv[1]);
	ImageType::Pointer dstImg = imageIO.ReadImageT(argv[2]);

	Timer timer;
	timer.start();
	computeMSE(srcImg, dstImg);
	timer.stop();

	cout << "Elapsed Time: " << timer.getElapsedTimeInMilliSec() << endl;

	timer.start();
	TransformType::Pointer transform = TransformType::New();
	transform->Rotate(0, 1, angle2rad(5));
	ImageType::Pointer movingImg1 = resampleImageITK(srcImg, transform);
	timer.stop();
	cout << "Elapsed Time: " << timer.getElapsedTimeInMilliSec() << endl;

	timer.start();
	ImageType::Pointer movingImg2 = resampleImageRaw(srcImg, srcImg, transform);
	timer.stop();
	cout << "Elapsed Time: " << timer.getElapsedTimeInMilliSec() << endl;
}

int main__(int argc, char* argv[]) {
	itkcmds::itkImageIO<ImageType> imageIO;
	ImageType::Pointer srcImg = imageIO.ReadImageT(argv[1]);
	ImageType::Pointer dstImg = imageIO.ReadImageT(argv[2]);

	TransformType::Pointer transform = TransformType::New();
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	MetricType::Pointer metric = MetricType::New();
	metric->SetFixedImage(dstImg);
	metric->SetFixedImageRegion(dstImg->GetLargestPossibleRegion());
	metric->SetMovingImage(srcImg);
	metric->SetInterpolator(interpolator);
	metric->SetUseAllPixels(false);
	metric->SetTransform(static_cast<MetricType::TransformPointer>(transform));
	metric->SetNumberOfThreads(8);
	metric->SetFixedImageSamplesIntensityThreshold(0.5);

	try {
		metric->Initialize();
	}
	catch (itkException &ex) {
		cout << ex << endl;
	}

	transform->Rotate(0, 1, 8 * itk::Math::pi / 180.);

	MetricType::ParametersType params = transform->GetParameters();
	MetricType::DerivativeType derivOut;

	double m = metric->GetValue(params);
	metric->GetDerivative(params, derivOut);

	cout << "Number of Pixels: " << metric->GetNumberOfPixelsCounted() << endl;
	cout << "Initial Params: " << params << endl;
	cout << "Metric Value: " << m << endl;
	cout << "Metric Derivative: " << derivOut << endl;

	return 0;
}

int run(int argc, char* argv[]) {

	MetricType::Pointer metric = MetricType::New();
	TransformType::Pointer transform = TransformType::New();
	OptimizerType::Pointer optimizer = OptimizerType::New();
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	RegistrationType::Pointer registration = RegistrationType::New();

	metric->SetUseAllPixels(true);

	optimizer->SetDebug(true);
	optimizer->SetGlobalWarningDisplay(true);

	registration->SetGlobalWarningDisplay(true);
	registration->SetDebug(true);
	registration->SetMetric(metric);
	registration->SetOptimizer(optimizer);
	registration->SetTransform(
			static_cast<RegistrationType::TransformPointer>(transform));
	registration->SetInterpolator(interpolator);

	itkcmds::itkImageIO<ImageType> imageIO;
	ImageType::Pointer srcImg = imageIO.ReadImageT(argv[1]);
	ImageType::Pointer dstImg = imageIO.ReadImageT(argv[2]);

	registration->SetFixedImage(dstImg);
	registration->SetMovingImage(srcImg);
	registration->SetFixedImageRegion(dstImg->GetLargestPossibleRegion());

	transform->Rotate(0, 1, 5 * itk::Math::pi / 180.);
	TransformType::ParametersType txParam = transform->GetParameters();

	RegistrationType::ParametersType initialParams(txParam);
	registration->SetInitialTransformParameters(txParam);

	optimizer->SetMaximumStepLength(4.0);
	optimizer->SetMinimumStepLength(0.01);
	optimizer->SetNumberOfIterations(200);

	try {
		registration->Update();
		std::cout << "Registration updated ..." << std::endl;
	}
	catch (itkException &ex) {
		std::cerr << ex;
		std::cerr << std::endl;
		return EXIT_FAILURE;
	}

	return 0;
}
