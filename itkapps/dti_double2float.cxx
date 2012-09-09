#include "itkImage.h"
#include "itkImageCommon.h"
#include <itkDiffusionTensor3D.h>
#include <itkCastImageFilter.h>
#include "itkImageRegionIteratorWithIndex.h"

#include "string"

using namespace std;

typedef itk::DiffusionTensor3D<double> DoubleTensorType ;
typedef itk::Image<DoubleTensorType, 3> DoubleImageType ;

typedef itk::DiffusionTensor3D<float> FloatTensorType ;
typedef itk::Image<FloatTensorType, 3> FloatImageType ;
typedef itk::ImageRegionIteratorWithIndex<DoubleImageType> DoubleImageIter;


int main(int argc, char* argv[]) {
	if (argc < 3) {
		cout << "usage: " << argv[0] << " in-dti out-dti" << endl;
		return 0;
	}

	int ret = 0;
	DoubleImageType::Pointer inImg = ReadImageT<DoubleImageType>(argv[1], ret);
	DoubleImageType::SizeType inSize = inImg->GetRequestedRegion().GetSize();

	FloatTensorType df0;
	df0.Fill(0.0);

	FloatImageType::Pointer outImg = NewImageT<FloatImageType>(inSize[0], inSize[1], inSize[2], df0);

	DoubleImageIter iter(inImg, inImg->GetRequestedRegion());
	for (; !iter.IsAtEnd(); ++iter) {
		DoubleImageType::IndexType idx = iter.GetIndex();

		DoubleTensorType dt = iter.Get();
		FloatTensorType df;

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				df(i,j) = dt(i,j);
			}
		}

		outImg->SetPixel(idx, df);
	}

	WriteImageT<FloatImageType>(argv[2], outImg);

}
