#include "itkImageCommon.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImage.h"
#include "iostream"

using namespace std;

typedef itk::Vector<double,3> VectorType;
typedef itk::Image<VectorType,3> VectorImageType;
typedef itk::Image<double,3> ImageType;
typedef itk::ImageRegionIteratorWithIndex<VectorImageType> IteratorType;

int main(int argc, char* argv[]) {
	if (argc < 6) {
		cout << "usage: " << argv[0] << " input-eigen-image input-cx input-cy input-cz output-h" << endl;
		return 0;
	}

	int ret = 0;
	ImageType::Pointer cx = ReadImageT<ImageType>(argv[2], ret);
	ImageType::Pointer cy = ReadImageT<ImageType>(argv[3], ret);
	ImageType::Pointer cz = ReadImageT<ImageType>(argv[4], ret);
	VectorImageType::Pointer eigenImg = ReadImageT<VectorImageType>(argv[1], ret);

	ImageType::Pointer hImg = NewImageT<ImageType>(cx);

	IteratorType iter(eigenImg, eigenImg->GetRequestedRegion());
	for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
		VectorImageType::IndexType idx = iter.GetIndex();
		VectorType vi = iter.Get();

		double vcx = cx->GetPixel(idx);
		double vcy = cy->GetPixel(idx);
		double vcz = cz->GetPixel(idx);

		double helicity = vi[0] * vcx + vi[1] * vcy + vi[2] * vcz;
		hImg->SetPixel(idx, helicity);
	}

	WriteImageT<ImageType>(argv[5], hImg);
}
