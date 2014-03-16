#include "itkImageCommon.h"

#include "iostream"
#include "itkImageRegionIteratorWithIndex.h"

using namespace std;

typedef itk::Image<double,3> ImageType;
typedef itk::Vector<double,3> VectorType;
typedef itk::Image<VectorType,3> VectorImageType;
typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;

int main(int argc, char* argv[]) {
  if (argc < 4) {
    cout << "usage: " << argv[0] << " input-phi-image input-theta-image output-vector-image" << endl;
    return 0;
  }

  int ret = 0;

  ImageType::Pointer phiImg = ReadImageT<ImageType>(argv[1], ret);
  ImageType::Pointer thetaImg  = ReadImageT<ImageType>(argv[2], ret);

  ImageType::RegionType phiRegion = phiImg->GetRequestedRegion();
  ImageType::SizeType phiSize = phiRegion.GetSize();

  VectorType zeroVector;
  zeroVector.Fill(0);

  VectorImageType::Pointer eigenImg = NewImageT<VectorImageType>(phiSize[0], phiSize[1], phiSize[2], zeroVector);

  IteratorType iter(phiImg, phiImg->GetRequestedRegion());
  for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
    ImageType::IndexType idx = iter.GetIndex();

    double phi = iter.Get();
    double theta = thetaImg->GetPixel(idx);
    double x = cos(phi) * cos(theta);
    double y = sin(phi) * cos(theta);
    double z = sin(theta);

    VectorType v;
    v[0] = x; v[1] = y; v[2] = z;

    eigenImg->SetPixel(idx, v);
  }

  WriteImageT<VectorImageType>(argv[3], eigenImg);
}
