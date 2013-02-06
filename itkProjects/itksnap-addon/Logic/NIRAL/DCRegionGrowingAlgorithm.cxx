#include "DCRegionGrowingAlgorithm.h"

// actual implementation of GrowingFunctionWithIndex
template <class LabelImageType>
bool niral::DCRegionGrowingAlgorithm<LabelImageType>::GrowingFunctionWithIndex(VectorImageType::IndexType i1, VectorImageType::PixelType v1,
      VectorImageType::IndexType i2, VectorImageType::PixelType v2, typename LabelImageType::PixelType &reason) {
  if (m_ROIimgLoaded) {
    unsigned short l = m_ROIimg->GetPixel(i2);
    if (l < 1) {
      return false;
    }
  }

  double m1 = mag3(v1);
  double m2 = mag3(v2);

  if (m1 == 0 || m2 == 0 || isnan(m1) || isnan(m2)) {
    return false;
  }

  double innerProduct = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
  double absIP = abs(innerProduct);

  bool included = absIP > m_DC;
  if (!included) {
    reason = 1;
    return false;
  }

  if (m_f1ImgLoaded) {
    double f1 = m_f1Img->GetPixel(i1);
    double f2 = m_f1Img->GetPixel(i2);
    if (m_f1Threshold[0] > f2 || f2 > m_f1Threshold[1]) {
			//cout << "F1: " << f2 << endl;
			included = false;
      reason = 2;
      return false;
		}
		if (abs(f1 - f2) > m_f1Delta) {
			included = false;
      reason = 3;
      return false;
		}
  }
  
  if (m_f2ImgLoaded) {
    double f1 = m_f2Img->GetPixel(i1);
    double f2 = m_f2Img->GetPixel(i2);
    if (m_f2Threshold[0] > f2 || f2 > m_f2Threshold[1]) {
			//cout << "F2: " << f2 << endl;
			included = false;
      reason = 4;
      return false;
		}
		if (abs(f1 - f2) > m_f2Delta) {
			included = false;
      reason = 5;
      return false;
		}
  }
  
  return included;
}

template <class LabelImageType>
void niral::DCRegionGrowingAlgorithm<LabelImageType>::PixelIncludedEvent(VectorImageType::IndexType nbrIdx) {
	++m_VoxelCount;
	if (m_VoxelCount % 10000 == 0) {
		cout << "Included Pixels: " << m_VoxelCount << endl;
	}
}

template <class LabelImageType>
void niral::DCRegionGrowingAlgorithm<LabelImageType>::InitializeEvent(IndexStackType seedStack) {
}
