#include "itkImageCommon.h"

#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkImageMomentsCalculator.h"
#include "itkNeighborhoodConnectedImageFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkNeighborhoodIterator.h"

using namespace std;

typedef itk::Image<unsigned int,3> ImageType;
typedef itk::LabelStatisticsImageFilter<ImageType,ImageType> StatFilter;
typedef itk::ConnectedComponentImageFilter<ImageType,ImageType> ConnFilter;
typedef itk::RelabelComponentImageFilter<ImageType,ImageType> RelabelFilter;
typedef itk::BinaryThresholdImageFilter<ImageType,ImageType> ThresholdFilter;
typedef itk::ImageMomentsCalculator<ImageType> MomentsCalculator;
typedef itk::NeighborhoodConnectedImageFilter<ImageType,ImageType> RegionGrowingFilter;
typedef itk::ImageRegionConstIteratorWithIndex<ImageType> ImageIterator;
typedef itk::NeighborhoodIterator<ImageType> NeighborIterator;

typedef itk::Vector<double,3> VectorType;
typedef itk::Matrix<double,3,3> MatrixType;

template<class TI1, class TI2=TI1, class TO=TI1>
class LabelMerge
{
	private:
		TI1 _label;

	public:
		LabelMerge() {};
		~LabelMerge() {};

		bool operator!=(const LabelMerge&) {
			return false;
		}
		bool operator==(const LabelMerge& other) {
			return !(*this != other);
		}
		void SetLabel(TI1 l) {
			_label = l;
		}
		inline TO operator()(const TI1& a, const TI2& b) const {
			if (b == 1) {
				return _label;
			}
			return a;
		}
};

ImageType::Pointer ExtractLargestSymmetricComponent(ImageType::Pointer im, int thLow, int thUp) {
	ThresholdFilter::Pointer tFilter = ThresholdFilter::New();
	tFilter->SetInput(im);
	tFilter->SetLowerThreshold(thLow);
	tFilter->SetUpperThreshold(thUp);
	tFilter->SetOutsideValue(0);
	tFilter->SetInsideValue(1);
	tFilter->Update();

	ImageType::Pointer thIm = tFilter->GetOutput();

	ConnFilter::Pointer conFilter = ConnFilter::New();
	conFilter->SetInput(thIm);
	conFilter->Update();

	RelabelFilter::Pointer relabelFilter = RelabelFilter::New();
	relabelFilter->SetInput(conFilter->GetOutput());
	relabelFilter->Update();

	ImageType::Pointer rlIm = relabelFilter->GetOutput();
	int nLabels = relabelFilter->GetNumberOfObjects();

	for (int i = 1; i <= nLabels; i++) {
		tFilter->SetInput(rlIm);
		tFilter->SetLowerThreshold(i);
		tFilter->SetUpperThreshold(i);
		tFilter->SetOutsideValue(0);
		tFilter->SetInsideValue(1);
		tFilter->Update();

		ImageType::Pointer ccIm = tFilter->GetOutput();
		MomentsCalculator::Pointer mc = MomentsCalculator::New();
		mc->SetImage(ccIm);
		mc->Compute();

		VectorType cog = mc->GetCenterOfGravity();
		cout << "Label: " << i << endl;
		cout << cog << endl;
	}
	return im;
}

ImageType::Pointer AccumulateConnectedComponents(ImageType::Pointer inIm, int nStart, int nEnd) {
	ImageType::SizeType szIm = inIm->GetRequestedRegion().GetSize();
	ImageType::Pointer outIm = NewImageT<ImageType>(szIm[0], szIm[1], szIm[2]);

	typedef itk::BinaryFunctorImageFilter<ImageType, ImageType, ImageType, LabelMerge<ImageType::PixelType> > MergeFilter;

	cout << "Max number of labels: " << nEnd << endl;
	for (int i = nStart; i < nEnd; i++) {
		ThresholdFilter::Pointer tFilter = ThresholdFilter::New();
		tFilter->SetInput(inIm);
		tFilter->SetLowerThreshold(i);
		tFilter->SetUpperThreshold(i);
		tFilter->SetInsideValue(1);
		tFilter->SetOutsideValue(0);
		tFilter->Update();

		ConnFilter::Pointer filter = ConnFilter::New();
		filter->SetInput(tFilter->GetOutput());
		filter->Update();

		RelabelFilter::Pointer relabelFilter = RelabelFilter::New();
		relabelFilter->SetInput(filter->GetOutput());
		relabelFilter->Update();

		MergeFilter::Pointer mFilter = MergeFilter::New();
		mFilter->GetFunctor().SetLabel((i+1));
		mFilter->SetInput1(outIm);
		mFilter->SetInput2(relabelFilter->GetOutput());
		mFilter->Update();

		outIm = mFilter->GetOutput();
	}

	return outIm;
}

ImageType::Pointer RegionGrowing(ImageType::Pointer inIm, ImageType::Pointer maskIm) {
	RegionGrowingFilter::Pointer filterPtr = RegionGrowingFilter::New();
	filterPtr->SetInput(inIm);
	filterPtr->SetLower(25);
	filterPtr->SetUpper(35);

	ImageIterator iter(maskIm, inIm->GetRequestedRegion());	
	int nSeeds = 0;
	while (!iter.IsAtEnd()) {
		if (iter.Get() == 1) {
			ImageType::IndexType idx = iter.GetIndex();
			ImageType::SizeType nbrRadius;
			nbrRadius.Fill(1);

			NeighborIterator nIter(nbrRadius, inIm, maskIm->GetRequestedRegion());
			nIter.SetLocation(idx);
			for (int j = 0; j < nIter.Size(); j++) {
				cout << nIter.GetPixel(j) << endl;
			}
			nSeeds ++;
		}
		++iter;
	}

	cout << "number of seeds: " << nSeeds << endl;

	ImageType::SizeType radius;
	radius[0] = radius[1] = radius[2] = 1;
	filterPtr->SetRadius(radius);
	filterPtr->SetReplaceValue(1000);
	filterPtr->Update();
	return filterPtr->GetOutput();
}


int main(int argc, char* argv[]) {
	if (argc < 3) {
		cout << "usage: " << argv[0] << " input-pcnn-image mask-image output-image" << endl;
		return 0;
	}

	int ret = 0;
	ImageType::Pointer inIm = ReadImageT<ImageType>(argv[1], ret);
	ImageType::Pointer maskIm = ReadImageT<ImageType>(argv[2], ret);

	StatFilter::Pointer stFilt = StatFilter::New();
	stFilt->SetInput(inIm);
	stFilt->SetLabelInput(inIm);
	stFilt->Update();
	int nLabels = stFilt->GetNumberOfLabels();

	ImageType::Pointer outIm;

	//outIm = AccumulateConnectedComponents(inIm, 1, nLabels);
	//outIm = ExtractLargestSymmetricComponent(outIm, 1, nLabels);	
	outIm = RegionGrowing(inIm, maskIm);
	WriteImageT<ImageType>(argv[3], outIm);

	return 0;
}
