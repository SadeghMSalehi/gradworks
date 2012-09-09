#include "itkImageCommon.h"
#include "itkBoxMeanImageFilter.h"
#include "iostream"

using namespace std;

typedef itk::Image<double,3> ImageType;
typedef itk::BoxMeanImageFilter<ImageType,ImageType> BoxMeanFilterType;

int main(int argc, char* argv[]) {
  if (argc < 3) {
    cout << "usage: " << argv[0] << " input-image output-image" << endl;
    return 0;
  }

  int ret = 0;
  ImageType::Pointer srcImg = ReadImageT<ImageType>(argv[1], ret);

  BoxMeanFilterType::Pointer filter = BoxMeanFilterType::New();
  filter->SetInput(srcImg);
  filter->Update();
  
  WriteImageT<ImageType>(argv[2], filter->GetOutput());
  return 0;
}
