#include "itkImageCommon.h"
#include "itkImageMomentsCalculator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "vtkMath.h"

#include "iostream"

using namespace std;

typedef itk::Image<unsigned short,3> ImageType;
typedef itk::ImageMomentsCalculator<ImageType> MomentsCalculator;
typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;

typedef itk::Vector<double,3> VectorType;
typedef itk::Matrix<double,3,3> MatrixType;

VectorType ComputePlaneNormal(MatrixType secondMoments) {

  double** sm, **evec;
  double *eval;

  sm = new double*[3];
  evec = new double*[3];
  eval = new double[3];


  for (int i = 0; i < 3; i++) {
    sm[i] = new double[3];
    evec[i] = new double[3];

    for (int j = 0; j < 3; j++) {
      sm[i][j] = secondMoments(i,j);
    }
  }

  vtkMath::Jacobi(sm, eval, evec);

  VectorType v;
  v[0] = evec[0][0];
  v[1] = evec[0][1];
  v[2] = evec[0][2];

  return v;
}


int main(int argc, char* argv[]) {
  if (argc < 6) {
    cout << "usage: " << argv[0] << " input-image output-image label-to-split old-label-splitted new-label-splitted " << endl;
    return 0;
  }

  int ret;
  ImageType::Pointer srcImg = ReadImageT<ImageType>(argv[1], ret);
	ImageType::SpacingType srcSpacing = srcImg->GetSpacing();

  MomentsCalculator::Pointer calc = MomentsCalculator::New();
  calc->SetImage(srcImg);
  calc->Compute();
  MatrixType secondMoments = calc->GetSecondMoments();

  VectorType centerOfGravity = calc->GetCenterOfGravity();
  VectorType planeNormal = ComputePlaneNormal(secondMoments);

  cout << "Plane normal: " << planeNormal << endl;
  cout << "Center of gravity: " << centerOfGravity << endl;

  IteratorType iter(srcImg, srcImg->GetRequestedRegion());

  int labelToSplit = atoi(argv[3]);
  int oldLabel = atoi(argv[4]);
  int newLabel = atoi(argv[5]);
  for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
    if (iter.Get() == labelToSplit) {
      ImageType::IndexType idx = iter.GetIndex();
      double v = planeNormal[0] * (srcSpacing[0]*idx[0] - centerOfGravity[0]) + planeNormal[1] * (srcSpacing[1]*idx[1] - centerOfGravity[1]) + planeNormal[2] * (srcSpacing[2]*idx[2] - centerOfGravity[2]);
      if (v > 0) {
        iter.Set(newLabel);
      } else {
        iter.Set(oldLabel);
			}
    }
  }

  WriteImageT<ImageType>(argv[2], srcImg);

}
