#include "itkImageCommon.h"
#include "PCNNSegmentationCLP.h"
#include "itkRescaleIntensityImageFilter.h"

typedef itk::Image<double,3> ImageType;

int main(int argc, char* argv[]) {
	PARSE_ARGS;

	int ret = 0;
	ImageType::Pointer imIn = ReadImageT<ImageType>(inputImage.c_str(), ret);


		
	return 0;
}
