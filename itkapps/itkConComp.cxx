#include "itkImageCommon.h"

#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkConnectedComponentFunctorImageFilter.h"

using namespace std;

typedef itk::Image<unsigned int,3> ImageType;
typedef itk::LabelStatisticsImageFilter<ImageType,ImageType> StatFilter;
typedef itk::ConnectedComponentImageFilter<ImageType,ImageType> ConnFilter;
typedef itk::RelabelComponentImageFilter<ImageType,ImageType> RelabelFilter;

template <class TInput>
class FixedIntensityPixelsFunctor {
	public:
		FixedIntensityPixelsFunctor() {};
		~FixedIntensityPixelsFunctor() {};

		bool operator!=(const FixedIntensityPixelsFunctor& other) {
			return (m_Intensity != other.m_Intensity);
		}

		bool operator==(const FixedIntensityPixelsFunctor& other) {
			return ((*this) != other);
		}

		bool operator()(const TInput& a, const TInput& b) {
			return ((a == m_Intensity) && (b == m_Intensity));
		}

		void SetIntensity(TInput intensity) {
			m_Intensity = intensity;
		}

		TInput GetIntensity() {
			return m_Intensity;
		}

	private:
		TInput m_Intensity = 0;
};

typedef itk::ConnectedComponentFunctorImageFilter<ImageType,ImageType,FixedIntensityPixelsFunctor<ImageType::PixelType>,ImageType> AllConnFilter;

int main(int argc, char* argv[]) {
	if (argc < 3) {
		cout << "usage: " << argv[0] << " input-pcnn-image output-image" << endl;
		return 0;
	}

	int ret = 0;
	ImageType::Pointer inIm = ReadImageT<ImageType>(argv[1], ret);

	StatFilter::Pointer stFilt = StatFilter::New();
	stFilt->SetInput(inIm);
	stFilt->SetLabelInput(inIm);
	stFilt->Update();

	ImageType::SizeType szIm = inIm->GetRequestedRegion().GetSize();
	ImageType::Pointer inOut = NewImageT<ImageType>(szIm[0], szIm[1], szIm[2]);

	int nLabels = stFilt->GetNumberOfLabels();

	for (int i = 1; i <= nLabels; i++) {
		AllConnFilter::Pointer filter = AllConnFilter::New();
		filter->GetFunctor().SetIntensity(i);
		filter->SetInput(inIm);
		filter->Update();

		RelabelFilter::Pointer relabelFilter = RelabelFilter::New();
		relabelFilter->SetInput(filter->GetOutput());
		relabelFilter->Update();
	}
}
