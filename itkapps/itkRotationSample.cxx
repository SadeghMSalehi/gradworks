#include "itkImageCommon.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "iostream"

using namespace std;

typedef itk::Vector<double,3> VectorType;
typedef itk::Image<VectorType,3> VectorImageType;
typedef itk::ImageRegionIteratorWithIndex<VectorImageType> IteratorType;

int main(int argc, char* argv[]) {
	VectorType v;
	v.Fill(0);

	VectorImageType::Pointer srcImg = NewImageT<VectorImageType>(100, 100, 100, v);
	IteratorType iter(srcImg, srcImg->GetRequestedRegion());
	
	for (; !iter.IsAtEnd(); ++iter) {
		VectorImageType::IndexType idx = iter.GetIndex();
		double x = idx[0] - 50;
		double y = idx[1] - 50;
		double z = idx[2] - 50;

		double r2 = x*x + y*y + z*z;

		if (r2 < 10*10) {
			v[0] = 0;
			v[1] = 0;
			v[2] = 0;
		} else if (r2 < 30*30) {
			v[0] = -y;
			v[1] = x;
			v[2] = 0;
		} else {
			v.Fill(0);
		}

		srcImg->SetPixel(idx, v);
	}

	WriteImageT<VectorImageType>(argv[1], srcImg);
}
