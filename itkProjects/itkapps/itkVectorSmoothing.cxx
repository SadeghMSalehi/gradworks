#include "itkImageCommon.h"
#include "itkImageRegionIteratorWithIndex.h" 
#include "itkPolarCoordinateFilter.h"
#include "itkMeanOnVectorImageFilter.h"
#include "math.h"

#include "iostream"

using namespace std;

typedef itk::ImageRegionIteratorWithIndex<VectorImageType> VectorIteratorType;
typedef itk::Image<double,3> ImageType;
typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;
typedef itk::CartesianToHalfSphericalCoordFilter SphericalCoordConversionFilterType;
typedef itk::SphericalToCartesianCoordFilter CartesianCoordConversionFilterType;
typedef itk::MeanOnVectorImageFilter<SphericalImageType,SphericalImageType> MeanFilterType;

int main(int argc, char* argv[]) {
  if (argc < 3) {
    cout << "usage: " << argv[0] << " input-vector-image output-vector-image" << endl;
    return 0;
  }

  int ret = 0;
  VectorImageType::Pointer srcImg = ReadImageT<VectorImageType>(argv[1], ret);

	VectorImageType::RegionType srcRegion = srcImg->GetRequestedRegion();
  VectorImageType::SizeType srcSize = srcImg->GetRequestedRegion().GetSize();

	cout << "Converting to spherical coord ..." << endl;
	SphericalCoordConversionFilterType::Pointer convFilter = SphericalCoordConversionFilterType::New();
	convFilter->SetInput(srcImg);
	convFilter->Update();
	SphericalImageType::Pointer sprImg = convFilter->GetOutput();
	WriteImageT<SphericalImageType>("sphericalImage.nrrd", sprImg);

	cout << "Computing mean of angles ..." << endl;
	MeanFilterType::Pointer meanFilter = MeanFilterType::New();
	meanFilter->SetInput(sprImg);
	SphericalImageType::Pointer sprMeanImg = meanFilter->GetOutput();
	WriteImageT<SphericalImageType>("sphericalMeanImage.nrrd", sprMeanImg);

	cout << "Converting back to cartesian coord ..." << endl;
	CartesianCoordConversionFilterType::Pointer convFilter2 = CartesianCoordConversionFilterType::New();
	convFilter2->SetInput(sprMeanImg);
	convFilter2->Update();
	WriteImageT<VectorImageType>(argv[2], convFilter2->GetOutput());

  return 0;
}
