#include "itkLabelMapMorphologyCLP.h"
#include "itkImageCommon.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryBallStructuringElement.h"

#include "iostream"
#include "string"

typedef itk::Image<unsigned short, 3> ImageType;
typedef itk::BinaryBallStructuringElement<unsigned int,3> KernelType;
typedef itk::BinaryDilateImageFilter<ImageType,ImageType,KernelType> DilationFilter;
typedef itk::BinaryErodeImageFilter<ImageType,ImageType,KernelType> ErosionFilter;


ImageType::Pointer RunDilation(ImageType::Pointer imIn, KernelType& ball, KernelType::PixelType dilateValue) {
	DilationFilter::Pointer df = DilationFilter::New();
	df->SetInput(imIn);
	df->SetKernel(ball);
	df->SetDilateValue(dilateValue);
	df->Update();
	return df->GetOutput();
}

ImageType::Pointer RunErosion(ImageType::Pointer imIn, KernelType& ball, KernelType::PixelType erosionValue) {
	ErosionFilter::Pointer ef = ErosionFilter::New();
	ef->SetInput(imIn);
	ef->SetKernel(ball);
	ef->SetErodeValue(erosionValue);
	ef->Update();
	return ef->GetOutput();
}

int main(int argc, char* argv[]) {
	PARSE_ARGS;

	int ret = 0;
	ImageType::Pointer imIn = ReadImageT<ImageType>(inputImage.c_str(), ret);

	KernelType ball;
	ball.SetRadius(size);
	ball.CreateStructuringElement();

	ImageType::Pointer temp;

	if (oper == "dilate") {
		temp = RunDilation(imIn, ball, 1);
	} else if (oper == "erode") {
		temp = RunErosion(imIn, ball, 1);
	} else if (oper == "close") {
		temp = RunDilation(imIn, ball, 1);
		temp = RunErosion(temp, ball, 1);
	} else if (oper == "open") {
		temp = RunErosion(imIn, ball, 1);
		temp = RunDilation(temp, ball, 1);
	}

	WriteImageT<ImageType>(outputImage.c_str(), temp);
}
