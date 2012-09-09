#include "itkImageCommon.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "iostream"

using namespace std;

typedef itk::Vector<double,3> VectorType;
typedef itk::Image<VectorType,3> VectorImageType;
typedef itk::Matrix<double,3,3> JacobianType;
typedef itk::Image<JacobianType, 3> JacobianImageType;
typedef itk::ImageRegionIteratorWithIndex<JacobianImageType> IteratorType;

int main(int argc, char* argv[]) {
  if (argc < 3) {
    cout << "usage: " << argv[0] << " input-vector-image output-jacobian-image" << endl;
    return 0;
  }

  int ret = 0;
  VectorImageType::Pointer srcImg = ReadImageT<VectorImageType>(argv[1], ret);
  VectorImageType::RegionType srcRegion = srcImg->GetRequestedRegion();
  VectorImageType::SizeType srcSize = srcRegion.GetSize();

  JacobianType zero;
  zero.Fill(0);
  JacobianImageType::Pointer outImg = NewImageT<JacobianImageType>(srcSize[0], srcSize[1], srcSize[2], zero);

  IteratorType iter(outImg, outImg->GetRequestedRegion());
  for (; !iter.IsAtEnd(); ++iter) {
    JacobianImageType::IndexType idx = iter.GetIndex();

    int x0 = idx[0] - 1;
    int x1 = idx[0] + 1;

    int y0 = idx[1] - 1;
    int y1 = idx[1] + 1;

    int z0 = idx[2] - 1;
    int z1 = idx[2] + 1;

    if (x0 < 0 || y0 < 0 || z0 < 0) {
      continue;
    }
    if (x1 >= (int) srcSize[0] || y1 >= (int) srcSize[1] || z1 >= (int) srcSize[2]) {
      continue;
    }
    
    VectorImageType::IndexType x0Idx = idx;
    VectorImageType::IndexType x1Idx = idx;

    x0Idx[0] = x0;
    x1Idx[0] = x1;

    VectorImageType::IndexType y0Idx = idx;
    VectorImageType::IndexType y1Idx = idx;

    y0Idx[1] = y0;
    y1Idx[1] = y1;

    VectorImageType::IndexType z0Idx = idx;
    VectorImageType::IndexType z1Idx = idx;

    z0Idx[2] = z0;
    z1Idx[2] = z1;

    VectorType x0v = srcImg->GetPixel(x0Idx);
    VectorType x1v = srcImg->GetPixel(x1Idx);
    VectorType y0v = srcImg->GetPixel(y0Idx);
    VectorType y1v = srcImg->GetPixel(y1Idx);
    VectorType z0v = srcImg->GetPixel(z0Idx);
    VectorType z1v = srcImg->GetPixel(z1Idx);

    VectorType dx = (x1v - x0v) / 2.0;
    VectorType dy = (y1v - y0v) / 2.0;
    VectorType dz = (z1v - z0v) / 2.0;

    JacobianType j;
    j(0,0) = dx[0];
    j(1,0) = dx[1];
    j(2,0) = dx[2];

    j(0,1) = dy[0];
    j(1,1) = dy[1];
    j(2,1) = dy[2];

    j(0,2) = dz[0];
    j(1,2) = dz[1];
    j(2,2) = dz[2];

    iter.Value() = j;
  }

  WriteImageT<JacobianImageType>(argv[2], outImg);
}
