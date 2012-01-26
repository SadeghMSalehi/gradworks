#ifndef __itkRegionGrowingAlgorithm_H__
#define __itkRegionGrowingAlgorithm_H__
#include "itkImageCommon.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "deque"
#include "limits.h"

using namespace std;

namespace niral {
  template <class T, class S>
  class RegionGrowingAlgorithm {
    public:
      typedef unsigned int uint32;
      
      typedef itk::Image<uint32, 3> VisitingImageType;
      typedef deque<typename T::IndexType> IndexStackType;


    private:
			int _numberOfIterations;
      IndexStackType _seedStack;

			// index to find the region of seed
			// record top-left and bottom-right cornder
			typename T::IndexType _minIndex;
			typename T::IndexType _maxIndex;
			
      typename T::RegionType _region;
      typename T::SizeType _size;

			// image to record region growing
			// out image is actual region growing results
			// propagation image will record the stage of region growing
      typename T::Pointer _inputImg;
      typename S::Pointer _seedImg;
      typename S::Pointer _outImg;
      typename S::Pointer _reasonForTerminationImage;
      typename S::Pointer _propagationImage;

			// to record the stage of region growing
      typename S::PixelType _curMarkingValue;
      
			// which pixels should visit or are visited?
			// recorded as a bitmap
      VisitingImageType::Pointer _visitImg;

      int _maxNumberOfPixels;
      int _numberOfPixels;

      RegionGrowingAlgorithm() {};

			// search seeds to begin region growing
      void findSeed() {
        _numberOfPixels = 0;

        typedef itk::ImageRegionConstIteratorWithIndex<S> IteratorType;
        typename S::SizeType seedSize;
        
        for (int j = 0; j < 3; j++) {
					if (_minIndex[j] > _maxIndex[j]) {
						seedSize[j] = 0;
					} else {
						seedSize[j] = _maxIndex[j] - _minIndex[j] + 1;
					}
				}
				
        typename S::RegionType seedRegion(_minIndex, seedSize);
        if (_numberOfIterations == 0) {
						seedRegion = _region;
				}
        
        for (IteratorType seedIter(_seedImg, seedRegion); !seedIter.IsAtEnd(); ++seedIter) {
          if (seedIter.Get() > 0) {
            typename S::IndexType seedIdx = seedIter.GetIndex();
            _seedStack.push_back(seedIdx);

            _outImg->SetPixel(seedIdx, seedIter.Get());
            _numberOfPixels ++;

						typename S::PixelType homeOrder = _propagationImage->GetPixel(seedIdx);
						if (homeOrder == 0) {
							homeOrder = static_cast<typename S::PixelType>(_numberOfPixels);
							_propagationImage->SetPixel(seedIdx, homeOrder);
						}
            // updateSeedExtent(seedIdx);
          }
        }
      }
      
			// update the extent of seed
      void updateSeedExtent(typename T::IndexType idx) {
				for (int i = 0; i < 3; i++) {
					_minIndex[i] = getMinValue(_minIndex[i], idx[i]);
					_maxIndex[i] = getMaxValue(_maxIndex[i], idx[i]);
				}
			}
			
      unsigned int getMaxValue(unsigned int x, unsigned int y) {
				return x > y ? x : y;
			}
			
			unsigned int getMinValue(unsigned int x, unsigned int y) {
				return x > y ? y : x;
			}


    public:
			// initialize variables
      RegionGrowingAlgorithm(typename T::Pointer inputImg, typename S::Pointer seedImg) {
        _curMarkingValue = 1;
        _inputImg = inputImg;
        _seedImg = seedImg;

        _region = _inputImg->GetLargestPossibleRegion();
        _size = _region.GetSize();
        
        // in case no region growing is processed
        // these index must designate the actual image size
        _maxIndex[0] = _maxIndex[1] = _maxIndex[2] = 0;
        _minIndex[0] = _size[0] - 1; 
        _minIndex[1] = _size[1] - 1; 
        _minIndex[2] = _size[2] - 1;

        _outImg = NewImageT<S>(_size[0], _size[1], _size[2]);
        CopyHeaderT<S>(_seedImg, _outImg);
        
        _visitImg = NewImageT<VisitingImageType>(_size[0], _size[1], _size[2]);
        _maxNumberOfPixels = INT_MAX;
        
        _reasonForTerminationImage = NewImageT<S>(_seedImg);
        _reasonForTerminationImage->FillBuffer(0);
        
        _propagationImage = NewImageT<S>(_seedImg);
        _propagationImage->FillBuffer(0);

        _numberOfPixels = 0;
        _numberOfIterations = 0;
      };

      virtual ~RegionGrowingAlgorithm() {
      };

      void SetMaxNumberOfPixels(int n) {
        _maxNumberOfPixels = n;
      }

			// utility function for visiting map
      static int offsetToIndex(typename T::OffsetType offset) {
        uint32 i = offset[0];
        uint32 j = offset[1];
        uint32 k = offset[2];

        return (i+1) + (j+1)*3 + (k+1)*9;
      }

      static typename T::OffsetType indexToOffset(int idx) {
        typename T::OffsetType offset;
        offset[2] = idx / 9 - 1;
        offset[1] = (idx % 9) / 3 - 1;
        offset[0] = idx % 3 - 1;
        return offset;
      }

      static void findNeighbors(typename T::Pointer inputImage, typename T::IndexType seed, typename T::PixelType* values) {
        for (int k = -1; k < 2; k+=1) {
          for (int j = -1; j < 2; j+=1) {
            for (int i = -1; i < 2; i+=1) {
              typename T::OffsetType offset;
              offset[0] = i; offset[1] = j; offset[2] = k;
              int idx = offsetToIndex(offset);
              typename T::IndexType nbrIdx = seed + offset;
              values[idx] = inputImage->GetPixel(nbrIdx);
            }
          }
        }
      }

      virtual void Init() {
        findSeed();
        InitializeEvent(_seedStack);
      }

      void SetOutputLabel(typename S::PixelType l) {
        _curMarkingValue = l;
      }

      /**
       * if already visited, return 0
       * if newly visiting, mark as visited and return 1
       */
      bool connect(VisitingImageType::Pointer visitImg, typename T::IndexType homeIdx, typename T::IndexType visitIdx) {
        // compare from home to visit
        uint32 visitorMark = visitImg->GetPixel(visitIdx);
        typename T::OffsetType offset = visitIdx - homeIdx;

        uint32 idx = offsetToIndex(offset);
        if ((visitorMark & (1 << idx)) > 0) {
          return false;
        }
        visitImg->SetPixel(visitIdx, visitorMark |= (1 << idx));

        // compare from visit to home
        uint32 homeMark = visitImg->GetPixel(homeIdx);
        typename T::OffsetType reverseOffset = homeIdx - visitIdx;
        uint32 reverseIdx = offsetToIndex(reverseOffset);
        /*
        if ((homeMark & (1 << reverseIdx)) > 0) {
          return false;
        }
        */
        visitImg->SetPixel(homeIdx, homeMark |= (1 << reverseIdx));
        return true;
      }

      virtual void DoRegionGrowing() {
        for (typename T::IndexType seedIdx = _seedStack.front();
						!_seedStack.empty() && _numberOfPixels < _maxNumberOfPixels; seedIdx = _seedStack.front(), _seedStack.pop_front()) {
          typename T::PixelType seedValue = _inputImg->GetPixel(seedIdx);

					// self connection - necessary??
          connect(_visitImg, seedIdx, seedIdx);
					updateSeedExtent(seedIdx);
					
          for (int i = 0; i < 27; i++) {
            typename T::OffsetType offset = indexToOffset(i);
            typename T::IndexType nbrIdx = seedIdx + offset;

            if (!_region.IsInside(nbrIdx)) {
              continue;
            }

            bool nbrMarked = _outImg->GetPixel(nbrIdx) > 0;
            bool nbrNotVisited = connect(_visitImg, seedIdx, nbrIdx);

            if (!nbrMarked && nbrNotVisited) {
              typename T::PixelType nbrValue = _inputImg->GetPixel(nbrIdx);
              typename S::PixelType reasonForTermination = 0;

              bool isIncluded = GrowingFunctionWithIndex(seedIdx, seedValue, nbrIdx, nbrValue, reasonForTermination);
              if (isIncluded) {
                _numberOfPixels++;
                _seedStack.push_back(nbrIdx);
                _outImg->SetPixel(nbrIdx, _curMarkingValue);
                updateSeedExtent(nbrIdx);
                PixelIncludedEvent(nbrIdx);
                _reasonForTerminationImage->SetPixel(nbrIdx, 0);
                _propagationImage->SetPixel(nbrIdx, _numberOfPixels);
              } else {
                _reasonForTerminationImage->SetPixel(nbrIdx, reasonForTermination);
              }
            }
          }
        }
      }
          
      virtual bool GrowingFunctionWithIndex(
          typename T::IndexType homeIdx, typename T::PixelType homeValue, 
          typename T::IndexType nbrIdx, typename T::PixelType nbrValue, typename S::PixelType& reason) {
        reason = 0;
        return false;
      }

			// must keep min and max index representing segmented region by previous region growing process
			// true if there is a segmentation result by previous region growing
			// false if there is no voxel contained in outImg
      virtual bool Reinitialize() {
				if (_numberOfPixels == 0) {
					return false;
				}
				
        // copy _outImg to _seedImg
        // ::DumpImageT<S>(_outImg, _seedImg);
        _seedImg = _outImg;

        // clean up visit map
        // _outImg->FillBuffer(0);
        _visitImg->FillBuffer(0);
        
        _numberOfPixels = 0;
        _numberOfIterations ++;

        // Initialize
        Init();
        return true;
      };      

      virtual void PixelIncludedEvent(typename T::IndexType nbrIdx) {
      };

      virtual void InitializeEvent(IndexStackType seedStack) {
      };
     
      bool IsOverPixels() {
        return _numberOfPixels >= _maxNumberOfPixels;
      }

      typename S::Pointer GetOutput() {
        return _outImg;
      }

      typename S::Pointer GetReasonForTermination() {
        return _reasonForTerminationImage;
      }
      
      typename S::Pointer GetPropagationImage() {
				_propagationImage->SetRequestedRegion(_propagationImage->GetLargestPossibleRegion());
				_propagationImage->SetBufferedRegion(_propagationImage->GetLargestPossibleRegion());
				return _propagationImage;
			}

  };
};

#endif
