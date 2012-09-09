#include "itkImageCommon.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "math.h"
#include "itkBoxMeanImageFilter.h"

#define HALF_PI (1.570796326794897)


#include "iostream"

using namespace std;

typedef itk::Vector<double,3> VectorType;
typedef itk::Image<VectorType,3> VectorImageType;
typedef itk::ImageRegionIteratorWithIndex<VectorImageType> VectorIteratorType;
typedef itk::Image<double,3> ImageType;
typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;
typedef itk::BoxMeanImageFilter<ImageType,ImageType> BoxMeanFilterType;


int main(int argc, char* argv[]) {
  if (argc < 4) {
    cout << "usage: " << argv[0] << " input-image output-phi-image output-theta-image" << endl;
    return 0;
  }

  int ret = 0;
  VectorImageType::Pointer srcImg = ReadImageT<VectorImageType>(argv[1], ret);

  VectorImageType::IndexType roiIdx;
  roiIdx[0] = 55; roiIdx[1] = 167; roiIdx[2] = 57;
  VectorImageType::SizeType roiSize;
  roiSize[0] = 150; roiSize[1] = 150; roiSize[2] = 150;

  VectorImageType::RegionType srcRegion(roiIdx, roiSize);
  VectorImageType::SizeType srcSize = srcImg->GetRequestedRegion().GetSize();

  ImageType::Pointer phiImg   = NewImageT<ImageType>(srcSize[0], srcSize[1], srcSize[2]);
  ImageType::Pointer thetaImg = NewImageT<ImageType>(srcSize[0], srcSize[1], srcSize[2]);

  ImageType::Pointer angleVarImg   = NewImageT<ImageType>(srcSize[0], srcSize[1], srcSize[2]);

  VectorIteratorType iter(srcImg, srcRegion);
  for (; !iter.IsAtEnd(); ++iter) {
    VectorType v = iter.Get();

    double phi = 0;
    if (v[0] < 0) {
      phi = atan2(-v[2], -v[0]);
    } else {
      phi = atan2(v[2], v[0]);
    } 
    double theta = acos(v[1]);

    phiImg->SetPixel(iter.GetIndex(), abs(phi));
    thetaImg->SetPixel(iter.GetIndex(), abs(HALF_PI - theta));
  }

  WriteImageT<ImageType>(argv[2], phiImg);
  WriteImageT<ImageType>(argv[3], thetaImg);

  BoxMeanFilterType::Pointer phiMeanFilter = BoxMeanFilterType::New();
  phiMeanFilter->SetInput(phiImg);
  phiMeanFilter->Update();

  BoxMeanFilterType::Pointer thetaMeanFilter = BoxMeanFilterType::New();
  thetaMeanFilter->SetInput(thetaImg);
  thetaMeanFilter->Update();


  return 0;
}
