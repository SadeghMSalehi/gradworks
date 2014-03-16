#include "itkImageCommon.h"
#include "itkExtractImageFilter.h"

typedef itk::ExtractImageFilter<ivImage, ivImage2> ivExtractFilter;

#include "iostream"

using namespace std;

int main(int argc, char* argv[]) {
	int ret;
	cout << "usage: " << argv[0] << " in-vector-image sliceNo out-vector-image" << endl;

	int sliceNo = atoi(argv[2]);
	ivImage::Pointer image = ReadVectorImage(argv[1], ret);

	ivImage::SizeType size;
	size = image->GetLargestPossibleRegion().GetSize();
	size[1] = 0;

	ivImage::IndexType index;
	index[0] = index[2] = 90;
	index[1] = sliceNo;

	ivExtractFilter::InputImageRegionType extractRegion(index, size);

	ivExtractFilter::Pointer extractFilter = ivExtractFilter::New();
	extractFilter->SetInput(image);
	extractFilter->SetExtractionRegion(extractRegion);
	extractFilter->Update();

	WriteImageT<ivImage2>(argv[3], extractFilter->GetOutput());
}
