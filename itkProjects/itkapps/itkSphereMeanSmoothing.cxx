#include "itkImageCommon.h"
#include "itkImageRegionIteratorWithIndex.h" 
#include "itkPolarCoordinateFilter.h"
#include "itkMeanOnVectorImageFilter.h"
#include "math.h"

#include "iostream"

using namespace std;

typedef itk::Vector<double,3> VectorType;
typedef itk::Image<VectorType,3> VectorImageType;
typedef itk::ImageRegionIteratorWithIndex<VectorImageType> VectorIteratorType;
typedef itk::MeanOnVectorImageFilter<VectorImageType,VectorImageType> MeanFilterType;

int main(int argc, char* argv[]) {
  if (argc < 3) {
    cout << "usage: " << argv[0] << " input-vector-image output-vector-image" << endl;
    return 0;
  }

  int ret = 0;
  VectorImageType::Pointer srcImg = ReadImageT<VectorImageType>(argv[1], ret);
	VectorImageType::RegionType srcRegion = srcImg->GetRequestedRegion();
  VectorImageType::SizeType srcSize = srcImg->GetRequestedRegion().GetSize();

	MeanFilterType::Pointer meanFilter = MeanFilterType::New();
	meanFilter->SetInput(srcImg);
  meanFilter->Update();
	VectorImageType::Pointer meanImg = meanFilter->GetOutput();
	WriteImageT<VectorImageType>(argv[2], meanImg);

  return 0;
}
