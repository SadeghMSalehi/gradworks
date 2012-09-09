#include "itkImageCommon.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "iostream"

using namespace std;

typedef itk::Vector<double,2> SphericalCoordType;
typedef itk::Image<SphericalCoordType,3> SphericalImageType;
typedef itk::Image<double,3> ImageType;
typedef itk::ImageRegionIteratorWithIndex<SphericalImageType> IteratorType;

int main(int argc, char* argv[]) {
  if (argc < 3) {
    cout << "usage: " << argv[0] << " input-vector-image output-image component-number" << endl;
    return 0;
  }

  int componentNumber = atoi(argv[3]);
  int ret = 0;
  SphericalImageType::Pointer srcImg = ReadImageT<SphericalImageType>(argv[1], ret);
  SphericalImageType::RegionType srcRegion = srcImg->GetRequestedRegion();
  SphericalImageType::SizeType srcSize = srcRegion.GetSize();
  ImageType::Pointer outImg = NewImageT<ImageType>(srcSize[0], srcSize[1], srcSize[2], 0);

  IteratorType iter(srcImg, srcRegion);
  for (; !iter.IsAtEnd(); ++iter) {
    SphericalCoordType v = iter.Get();
    outImg->SetPixel(iter.GetIndex(), v[componentNumber]); 
  }

  WriteImageT<ImageType>(argv[2], outImg);

}
