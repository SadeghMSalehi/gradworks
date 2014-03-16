#include "itkImageCommon.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "iostream"
#include "string"
#include "vector"
#include "dtiCaseMeasureCLP.h"

using namespace std;

typedef itk::Image<unsigned short,3> ImageType;
typedef itk::Image<double,3> OutputImageType;
typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;

void writeImage(string f, OutputImageType::Pointer p) {
	if (f != "") {
		WriteImageT<OutputImageType>(f.c_str(), p);
	}
}

int main(int argc, char* argv[]) {
	PARSE_ARGS;

  int ret = 0;

  ImageType::Pointer l1 = ReadImageT<ImageType>(inputLambda1.c_str(), ret);
  ImageType::Pointer l2 = ReadImageT<ImageType>(inputLambda2.c_str(), ret);
  ImageType::Pointer l3 = ReadImageT<ImageType>(inputLambda3.c_str(), ret);

  ImageType::RegionType region = l1->GetRequestedRegion();
  ImageType::SizeType sz = region.GetSize();

  OutputImageType::Pointer clImg = NewImageT<OutputImageType>(sz[0], sz[1], sz[2]);
  OutputImageType::Pointer cpImg = NewImageT<OutputImageType>(sz[0], sz[1], sz[2]);
  OutputImageType::Pointer csImg = NewImageT<OutputImageType>(sz[0], sz[1], sz[2]);
  OutputImageType::Pointer raImg = NewImageT<OutputImageType>(sz[0], sz[1], sz[2]);
  OutputImageType::Pointer faImg = NewImageT<OutputImageType>(sz[0], sz[1], sz[2]);

  IteratorType iter(l1, l1->GetRequestedRegion());

  for (; !iter.IsAtEnd(); ++iter) {
    ImageType::IndexType idx = iter.GetIndex();
    ImageType::PixelType v1 = l1->GetPixel(idx);
    ImageType::PixelType v2 = l2->GetPixel(idx);
    ImageType::PixelType v3 = l3->GetPixel(idx);

    if (!(v1 >= v2 && v2 >= v3)) {
      cout << "Wrong order of lambdas: " << v1 << "," << v2 << "," << v3 << endl;
    }

    double vs = double(v1 + v2 + v3);
    double cl = (v1 - v2) / vs;
    double cp = 2.0 * (v2 - v3) / vs;
    double cs = 3.0 * v3 / vs;

    clImg->SetPixel(idx, cl);
    cpImg->SetPixel(idx, cp);
    csImg->SetPixel(idx, cs);

		if (outputFA != "") {
			double fa = 1 / sqrt(2) * sqrt((v1-v2)*(v1-v2) + (v2-v3)*(v2-v3) + (v3-v1)*(v3-v1)) / sqrt(v1*v1 + v2*v2 + v3*v3);
			faImg->SetPixel(idx, fa);
		}

		if (outputRA != "") {
			double ra = sqrt(3) * sqrt( (v1-vs/3)*(v1-vs/3) + (v2-vs/3)*(v2-vs/3) + (v3-vs/3)*(v3-vs/3) ) / vs;
			raImg->SetPixel(idx, ra);
		}
  }

	writeImage(outputCL, clImg);
	writeImage(outputCP, cpImg);
	writeImage(outputCS, csImg);
	writeImage(outputFA, faImg);
	writeImage(outputRA, raImg);
}
