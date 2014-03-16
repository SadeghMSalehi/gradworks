#include "itkImageCommon.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkGradientImageFilter.h"
#include "iostream"
#include "math.h"

using namespace std;

typedef itk::Vector<double,3> VectorType;
typedef itk::Image<VectorType,3> VectorImageType;
typedef itk::Image<double,3> ImageType;
typedef itk::ImageRegionIteratorWithIndex<ImageType> ImageIteratorType;

typedef itk::Matrix<double,3,3> JacobianType;
typedef itk::Image<JacobianType,3> JacobianImageType;
typedef itk::ImageRegionIterator<JacobianImageType> IteratorType;

typedef itk::Vector<VectorType,3> VectorVectorType;
typedef itk::Image<VectorVectorType,3> VectorVectorImageType;

typedef itk::GradientImageFilter<VectorImageType, VectorType, VectorVectorType> VectorGradientFilter;


int oldCode1(int argc, char* argv[]) { if (argc < 3) { cout << "usage: " << argv[0] << " input-jacobian output-curl-image" << endl;
    return 0;
  }

  int ret = 0;
  JacobianImageType::Pointer srcImg = ReadImageT<JacobianImageType>(argv[1], ret);
  IteratorType iter(srcImg, srcImg->GetRequestedRegion());

  JacobianImageType::RegionType srcRegion = srcImg->GetRequestedRegion();
  JacobianImageType::SizeType srcSize = srcRegion.GetSize();

  VectorType zeroVector;
  zeroVector.Fill(0.0);

  VectorImageType::Pointer outImg = NewImageT<VectorImageType>(srcSize[0], srcSize[1], srcSize[2], zeroVector);
  for (; !iter.IsAtEnd(); ++iter) {
    JacobianType t = iter.Get();
    double cx = t(1,2) - t(2,1);
    double cy = t(0,2) - t(2,0);
    double cz = t(0,1) - t(1,0);

    VectorType v;
    v[0] = cx; 
    v[1] = cy;
    v[2] = cz;

    outImg->SetPixel(iter.GetIndex(), v);
  }

  WriteImageT<VectorImageType>(argv[2], outImg);
  return 0;
}

int main(int argc, char* argv[]) {
	if (argc < 7) {
		cout << "usage: " << argv[0] << " input-dx input-dy input-dz output-curl-x output-curl-y output-curl-z" << endl;
		return 0;
	}

	int ret = 0;
	VectorImageType::Pointer dxImg = ReadImageT<VectorImageType>(argv[1], ret);
	VectorImageType::Pointer dyImg = ReadImageT<VectorImageType>(argv[2], ret);
	VectorImageType::Pointer dzImg = ReadImageT<VectorImageType>(argv[3], ret);

	VectorImageType::SizeType imgSize = dxImg->GetRequestedRegion().GetSize();

	ImageType::Pointer cxImg = NewImageT<ImageType>(imgSize[0], imgSize[1], imgSize[2]);
	ImageType::Pointer cyImg = NewImageT<ImageType>(imgSize[0], imgSize[1], imgSize[2]);
	ImageType::Pointer czImg = NewImageT<ImageType>(imgSize[0], imgSize[1], imgSize[2]);
	ImageType::Pointer cmImg = NewImageT<ImageType>(imgSize[0], imgSize[1], imgSize[2]);

	ImageIteratorType iter(cxImg, cxImg->GetRequestedRegion());
	for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
		ImageType::IndexType idx = iter.GetIndex();

		VectorType vx = dxImg->GetPixel(idx);
		VectorType vy = dyImg->GetPixel(idx);
		VectorType vz = dzImg->GetPixel(idx);

		double cx = vz[1] - vy[2];
		double cy = vx[2] - vz[0];
		double cz = vy[0] - vx[1];

		cxImg->SetPixel(idx, cx);
		cyImg->SetPixel(idx, cy);
		czImg->SetPixel(idx, cz);

		double cm = sqrt(cx*cx + cy*cy + cz*cz);
		cmImg->SetPixel(idx, cm);
	}

	WriteImageT<ImageType>(argv[4], cxImg);
	WriteImageT<ImageType>(argv[5], cyImg);
	WriteImageT<ImageType>(argv[6], czImg);
	WriteImageT<ImageType>("curlMagnitude.nrrd", cmImg);

	return 0;
}

/*
int main(int argc, char* argv[]) {
	if (argc < 3) {
		cout << "usage: " << argv[0] << " input-vector output-vector " << endl;
		return;
	}

	int ret = 0;
	VectorImageType::Pointer eigenImg = ReadImageT<VectorImageType>(argv[1], ret);
	VectorGradientFilter::Pointer gradFilter = VectorGradientFilter::New();
	gradFilter->SetInput(eigenImg);
	gradFilter->Update();
	VectorVectorImageType::Pointer jacobImg = gradFilter->GetOutput();


}
<<<<<<< itkCurlFilter.cxx

int main(int argc, char* argv[]) {
  if (argc < 3) {
    cout << "usage: " << argv[0] << " input-vector output-jacobian" << endl;
  }

	VectorImageType::Pointer srcImg = ReadImageT<VectorImageType>(argv[1], ret);

}
=======
*/
