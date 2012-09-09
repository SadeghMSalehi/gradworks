#include "iostream"
#include "stack"
#include "itkImageCommon.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkRegionGrowingAlgorithm.h"
#include "math.h"

#define PI (3.141592653589793)
#define mag3(x) (x[0]*x[0]+x[1]*x[1]+x[2]*x[2])

typedef itk::Vector<double,3> EigenVectorType;
typedef itk::Image<EigenVectorType,3> EigenImageType;
typedef itk::Image<double,3> ImageType;

ImageType::Pointer computeInnerProduct(EigenImageType::Pointer inputImg) {
  EigenImageType::RegionType region = inputImg->GetRequestedRegion();
  EigenImageType::SizeType size = region.GetSize();
  cout << "Input image size: " << size << endl;

  ImageType::Pointer outImg = NewImageT<ImageType>(size[0], size[1], size[2]);
  typedef itk::ImageRegionIteratorWithIndex<ImageType> IterType;
  IterType iter(outImg, region);

  int count = 0, previousProgress = 0;
  int nbrIdx[6] = { 4, 10, 12, 14, 16, 22 };

  for (; !iter.IsAtEnd(); ++iter) {
    ImageType::IndexType idx = iter.GetIndex();
    EigenImageType::PixelType curVal = inputImg->GetPixel(idx);

    double maxIP = 0;
    for (int j = 0; j < 6; j++) {
      int i = nbrIdx[j];

      ++count;
      int currentProgress = round((count * 100.0 / 27.0) / double(size[0] * size[1] * size[2]));
      if (currentProgress > previousProgress) {
        cout << currentProgress << "%" << endl;
        previousProgress = currentProgress;
      }

      if (i == 13) {
        continue;
      }

      ImageType::OffsetType offset = niral::RegionGrowingAlgorithm<ImageType,ImageType>::indexToOffset(i);
      ImageType::IndexType nbrIdx = idx + offset;

      if (!region.IsInside(nbrIdx)) {
        continue;
      }

      EigenImageType::PixelType nbrVal = inputImg->GetPixel(nbrIdx);
      if (isnan(mag3(curVal)) || isnan(mag3(nbrVal)) || mag3(curVal) == 0 || mag3(nbrVal) == 0) {
        continue;
      }

      double innerProduct = abs(curVal[0] * nbrVal[0] + curVal[1] * nbrVal[1] + curVal[2] * nbrVal[2]);
      if (innerProduct > maxIP) {
        maxIP = innerProduct;
      }
    }

    outImg->SetPixel(idx, acos(maxIP) * 180 / PI);
  }
  return outImg;
}

int main(int argc, char* argv[]) {
	if (argc < 3) {
		cout << "Usage: " << argv[0] << " input-image output-image" << endl;
		return 0;
	}

	int ret = 0;

	cout << "Reading " << argv[1];
	cout.flush();
	EigenImageType::Pointer inputImg = ReadImageT<EigenImageType>(argv[1], ret);
	cout << " done." << endl;

	WriteImageT<ImageType>(argv[2], computeInnerProduct(inputImg));
}
