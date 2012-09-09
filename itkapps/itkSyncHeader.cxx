#include "itkImageCommon.h"
#include "iostream"

using namespace std;

typedef itk::Image<unsigned short,3> ImageType;

int main(int argc, char* argv[]) {
  if (argc < 2) {
    cout << "usage: " << argv[0] << " src-img dst-img1 dst-img2 ... " << endl;
    return 0;
  }

  int ret = 0;
  ImageType::Pointer srcImg = ReadImageT<ImageType>(argv[1], ret);
  cout << "Source Spacing: " << srcImg->GetSpacing() << endl;
  cout << "Source Origin:  " << srcImg->GetOrigin() << endl;
  
  for (int i = 2; i < argc; i++) {
    ImageType::Pointer dstImg = ReadImageT<ImageType>(argv[2], ret);
    dstImg->SetSpacing(srcImg->GetSpacing());
    dstImg->SetOrigin(srcImg->GetOrigin());
    WriteImageT<ImageType>(argv[2], dstImg);
  }
}
