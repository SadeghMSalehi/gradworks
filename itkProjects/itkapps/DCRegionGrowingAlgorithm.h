#ifndef __DCRegionGrowingAlgorithm_H__
#define __DCRegionGrowingAlgorithm_H__

#include "itkOrientedImage.h"
#include "itkRegionGrowingAlgorithm.h"
#include "string"

namespace niral {
  typedef itk::Vector<double,3> VectorType;
  typedef itk::OrientedImage<VectorType,3> VectorImageType;
  typedef itk::OrientedImage<double,3> FAImageType;

  class GlobalStatistics {
    private:
      double _sumX;
      double _sumX2;
      int _n;
    
    public:
      GlobalStatistics() {
        _sumX = 0;
        _sumX2 = 0;
        _n = 0;
      }

      virtual ~GlobalStatistics() {}

      void Init() {
        _n = 0;
        _sumX = 0;
        _sumX2 = 0;
      }

      int GetN() { return _n; }
      
      void Add(double x) {
        _sumX += x;
        _sumX2 += (x*x);
        _n++;
      }

      double GetMean() {
        if (_n > 1) {
          return _sumX / double(_n);
        } else {
          return 0;
        }
      }

      double GetVar() {
        if (_n > 1) {
          double mean = GetMean();
          return (_sumX2 / double(_n) - mean * mean);
        } else {
          return 0;
        }
      }
      
      double GetStd() {
        double var = GetVar();
        return sqrt(var);
      }
      
      void Print() {
        cout << "N=" << _n << " Mean=" << GetMean() << " STD=" << GetStd() << endl;
      }
  };

	// specialized region growing algorithm for DTI
	// this class takes two feature volumes as inputs, feature1 image and feature2 image
	// then return its growing question with GrowingFunctionWithIndex overriden function
	// at each iteration (stage) of region growing, Reinitialize() must be called to clean up current region growing and update its seed
	template<class LabelImageType>
  class DCRegionGrowingAlgorithm : public RegionGrowingAlgorithm<VectorImageType,LabelImageType>
  {
    public:
    
			typedef typename RegionGrowingAlgorithm<VectorImageType,LabelImageType>::IndexStackType IndexStackType;
			
      DCRegionGrowingAlgorithm(VectorImageType::Pointer srcImg, typename LabelImageType::Pointer seedImg) 
          : RegionGrowingAlgorithm<VectorImageType,LabelImageType>(srcImg, seedImg) {
        m_f1Img = m_f2Img = NULL;
        m_ROIimg = NULL;
        m_DC = 2;
        m_InitialStd = 0;
        m_CutoffStd = 1;
        m_VoxelCount = 0;
        m_ROIimgLoaded = false;
        m_f1ImgLoaded = m_f2ImgLoaded = false;
        m_f1Threshold[0] = m_f1Threshold[1] = m_f1Delta = 0;
        m_f2Threshold[0] = m_f2Threshold[1] = m_f2Delta = 0;
      };

      virtual bool GrowingFunctionWithIndex(VectorImageType::IndexType, VectorImageType::PixelType, 
            VectorImageType::IndexType, VectorImageType::PixelType, typename LabelImageType::PixelType &reason);

      virtual void PixelIncludedEvent(VectorImageType::IndexType);
      virtual void InitializeEvent(IndexStackType);

      virtual bool Reinitialize() {
        bool keepGoing = RegionGrowingAlgorithm<VectorImageType,LabelImageType>::Reinitialize();
        m_Stat.Init();
        return keepGoing;
      }


      void SetROI(typename LabelImageType::Pointer roiImg) {
        m_ROIimg = roiImg;
        m_ROIimgLoaded = true;
      }

      void SetDC(double a) { m_DC = a; };
      void SetInitialStd(double s) { m_InitialStd = s; };
      void SetCutoffStd(double s) { m_CutoffStd = s; };
      
      void SetFeature1Image(FAImageType::Pointer fa) {
        m_f1Img = fa;
        m_f1ImgLoaded = true;
      }
      
      void SetFeature2Image(FAImageType::Pointer fa) {
        m_f2Img = fa;
        m_f2ImgLoaded = true;
      }      
    
      void SetFeature1Values(double min, double max, double delta) { 
				m_f1Threshold[0] = min;
				m_f1Threshold[1] = max;
				m_f1Delta = delta;
				cout << "Feature 1 threshold (" << m_f1Threshold[0] << ", " << m_f1Threshold[1] << ")" << endl;
			};
      
      void SetFeature2Values(double min, double max, double delta) { 
				m_f2Threshold[0] = min;
				m_f2Threshold[1] = max;
				m_f2Delta = delta;
				cout << "Feature 2 threshold (" << m_f2Threshold[0] << ", " << m_f2Threshold[1] << ")" << endl;
			};
			
      int GetVoxelCount() { return m_VoxelCount; };

    protected:
      inline double mag3(VectorType x) {
        return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
      }

    private:
      GlobalStatistics m_Stat;

      typename LabelImageType::Pointer m_ROIimg;

      double m_DC;
      double m_InitialStd;
      double m_CutoffStd;

      int m_VoxelCount;
      bool m_ROIimgLoaded;
      
      bool m_f1ImgLoaded;
      bool m_f2ImgLoaded;

      FAImageType::Pointer m_f1Img;
      FAImageType::Pointer m_f2Img;
      
      double m_f1Threshold[2];
      double m_f2Threshold[2];
      double m_f1Delta, m_f2Delta;
  };
};

#include "DCRegionGrowingAlgorithm.cxx"
#endif
