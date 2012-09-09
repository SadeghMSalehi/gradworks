#include "itkImageCommon.h"
#include "itkImageMomentsCalculator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "vtkMath.h"

#include "iostream"

using namespace std;

typedef itk::Image<double,3> ImageType;
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
  if (argc < 2) {
    cout << "usage: " << argv[0] << " input-image " << endl;
    return 0;
  }

  int ret;
  ImageType::Pointer srcImg = ReadImageT<ImageType>(argv[1], ret);
	ImageType::SpacingType srcSpacing = srcImg->GetSpacing();

  MomentsCalculator::Pointer calc = MomentsCalculator::New();
  calc->SetImage(srcImg);
  calc->Compute();
  MatrixType secondMoments = calc->GetSecondMoments();

	cout << secondMoments << endl;

  VectorType centerOfGravity = calc->GetCenterOfGravity();

  cout << "Center of gravity: " << centerOfGravity << endl;
}
