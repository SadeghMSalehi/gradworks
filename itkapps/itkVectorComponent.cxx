#include "itkImageCommon.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "iostream"

using namespace std;

int main(int argc, char* argv[]) {
  if (argc < 5) {
    cout << "usage: " << argv[0] << " input-vector-image output-component1 output-component2 output-component3 " << endl;
    return 0;
  }

  int numberOfComponents = 3;
	typedef itk::Vector<double,3> VectorType;
	typedef itk::Image<VectorType,3> VectorImageType;
	typedef itk::Image<double,3> ImageType;
	typedef itk::ImageRegionIteratorWithIndex<VectorImageType> IteratorType;

  int ret = 0;
  VectorImageType::Pointer srcImg = ReadImageT<VectorImageType>(argv[1], ret);
  VectorImageType::RegionType srcRegion = srcImg->GetRequestedRegion();
  VectorImageType::SizeType srcSize = srcRegion.GetSize();

	for (int i = 0; i < numberOfComponents; i++) {
		ImageType::Pointer outImg = NewImageT<ImageType>(srcSize[0], srcSize[1], srcSize[2], 0);

		IteratorType iter(srcImg, srcRegion);
		for (; !iter.IsAtEnd(); ++iter) {
			VectorType v = iter.Get();
			outImg->SetPixel(iter.GetIndex(), v[i]); 
		}

		WriteImageT<ImageType>(argv[i + 2], outImg);
	}

}
