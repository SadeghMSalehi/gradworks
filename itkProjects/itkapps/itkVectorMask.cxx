#include "itkImageCommon.h"
#include "itkImageRegionIterator.h"
#include "iostream"

using namespace std;

typedef itk::Vector<double,3> VectorType;
typedef itk::Image<VectorType,3> VectorImageType;
typedef itk::Image<unsigned short, 3> LabelImageType;
typedef itk::ImageRegionIterator<VectorImageType> VectorImageIteratorType;
typedef itk::ImageRegionIterator<LabelImageType> LabelImageIteratorType;

int main(int argc, char* argv[]) {
	int ret = 0;

	VectorImageType::Pointer dtivImage = ReadImageT<VectorImageType>(argv[1], ret);
	LabelImageType::Pointer seedImage = ReadImageT<LabelImageType>(argv[2], ret);

	VectorImageIteratorType it1(dtivImage, dtivImage->GetRequestedRegion());
	LabelImageIteratorType  it2(seedImage, seedImage->GetRequestedRegion());

	it1.GoToBegin();
	it2.GoToBegin();

	VectorType zeroVector;
  zeroVector.Fill(0);

	while (!it1.IsAtEnd()) {
		if (it2.Get() == 0) { 
				it1.Set(zeroVector);	
		} else {
      cout << it1.Get() << endl;
    }
		++it1;
		++it2;
	}

	WriteImageT<VectorImageType>(argv[3], dtivImage, false);
	return 0;
}
