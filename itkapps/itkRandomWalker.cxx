#include "itkImageCommon.h"
#include "itkImageSliceIteratorWithIndex.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "math.h"
#include "stdio.h"
#include "time.h"
#include "spMatrix.h"
#include "iostream"
#include "vector"


#include <vector>

using namespace std;

typedef itk::Image<double,3> ImageType;
typedef itk::Image<short,3> LabelMapType;
typedef itk::Image<unsigned int,3> OrderMapType;

typedef itk::ImageSliceIteratorWithIndex<OrderMapType> SliceIteratorType;

template <class TImage>
class ImageAccessor {
  public:
		ImageAccessor() {
		};

    ImageAccessor(typename TImage::Pointer img) {
      _img = img;
      _imgSz = img->GetBufferedRegion().GetSize();
    }

		void SetImage(typename TImage::Pointer img) {
			_img = img;
      _imgSz = img->GetBufferedRegion().GetSize();
		}

    int SliceIndex2Sub(typename TImage::IndexType idx) {
      return idx[0] + idx[2] * _imgSz[0];
    }

    typename TImage::IndexType Sub2SliceIndex(int slice, int sub) {
      typename TImage::IndexType idx;
      idx[0] = sub % _imgSz[0];
      idx[1] = slice;
      idx[2] = int(sub / _imgSz[0]);
      return idx;
    }
  private:
    typename TImage::Pointer _img;
    typename TImage::SizeType _imgSz;
};

template<class TImage>
class NeighborIterator {
  public:
    NeighborIterator(typename TImage::Pointer img) {
      _img = img;
      _curLoc = 0;
      _xyOffsets[0][0] = -1; _xyOffsets[0][1] = 0; _xyOffsets[0][2] = -1;
			_xyOffsets[1][0] = 0; _xyOffsets[1][1] = 0; _xyOffsets[1][2] = -1;
			_xyOffsets[2][0] = 1; _xyOffsets[2][1] = 0; _xyOffsets[2][2] = -1;
			_xyOffsets[3][0] = -1; _xyOffsets[3][1] = 0; _xyOffsets[3][2] = 0;
			_xyOffsets[4][0] = 1; _xyOffsets[4][1] = 0; _xyOffsets[4][2] = 0;
			_xyOffsets[5][0] = -1; _xyOffsets[5][1] = 0; _xyOffsets[5][2] = 1;
			_xyOffsets[6][0] = 0; _xyOffsets[6][1] = 0; _xyOffsets[6][2] = 1;
			_xyOffsets[7][0] = 1; _xyOffsets[7][1] = 0; _xyOffsets[7][2] = 1;
    }   

    void GoTo(typename TImage::IndexType idx) {
      _curIdx = idx;
      _curLoc = 0;
      ++(*this);
    }   
    bool IsAtEnd() {
      return _curLoc == 8;
    }
    void operator++() {
      do {
        _curLoc ++;
      } while (!_img->GetBufferedRegion().IsInside(GetIndex()) && _curLoc < 8);
    }
    typename TImage::IndexType GetIndex() {
      typename TImage::IndexType idx = _curIdx;
      for (int i = 0; i < 3; i++) {
        idx[i] = idx[i] + _xyOffsets[_curLoc][i];
      }
      return idx;
    }
    typename TImage::PixelType Get() {
      return _img->GetPixel(GetIndex());
    }
  private:
    typename TImage::Pointer _img;
    typename TImage::IndexType _curIdx;
    int _curLoc;
    int _xyOffsets[8][3];
};

class RandomWalker {
	private:
		ImageType::Pointer _srcImg;
		LabelMapType::Pointer _labelImg;
		OrderMapType::Pointer _orderImg;

		ImageType::Pointer _probImg;

		ImageType::SizeType _imgSz;
		ImageType::RegionType _imgRegion;

		spMatrix Lu, B;
		spError err;

		vector<ImageType::IndexType> _label0Index;
		vector<ImageType::IndexType> _label1Index;

		ImageAccessor<LabelMapType> _labelAccr;
		int _labelCount[255];

	public:
		RandomWalker(ImageType::Pointer srcImg, LabelMapType::Pointer labelImg) {
			_srcImg = srcImg;
			_labelImg = labelImg;

			_imgRegion = srcImg->GetBufferedRegion();
			_imgSz = _imgRegion.GetSize();

			_labelAccr.SetImage(labelImg);
			_probImg = NewImageT<ImageType>(_imgSz[0], _imgSz[1], _imgSz[2]);
			_orderImg = NewImageT<OrderMapType>(_imgSz[0], _imgSz[1], _imgSz[2], 0);

			Lu = spCreate(1000000, 0, &err);
			cout << "Lu created : " << err << endl;
			B = spCreate(1000000, 0, &err);
			cout << "B created : " << err << endl;
		}

		ImageType::Pointer GetProbabilityImage(LabelMapType::PixelType label) {
			return _probImg;
		}

		double ComputeWeight(ImageType::IndexType curIdx, ImageType::IndexType nbrIdx, ImageType::PixelType curVal, ImageType::PixelType nbrVal) {
				double weight = exp(-0.5*(curVal - nbrVal)*(curVal - nbrVal));
					return weight;
		}

		void Run() {
			SliceIteratorType iter(_orderImg, _imgRegion);
			iter.SetFirstDirection(0);
			iter.SetSecondDirection(2);
			time_t now;

			while (!iter.IsAtEnd()) {
				spClear(Lu);

				for (int i = 0; i < 255; i++) {
					_labelCount[i] = 0;
				}

				while (!iter.IsAtEndOfSlice()) {
					while (!iter.IsAtEndOfLine()) {
						OrderMapType::IndexType curIdx = iter.GetIndex();
						OrderMapType::PixelType curId = iter.Get();
						LabelMapType::PixelType curLabel = _labelImg->GetPixel(curIdx);
						_orderImg->SetPixel(curIdx, _labelCount[curLabel]);
						_labelCount[curLabel] ++;
						++iter;
					}
					iter.NextLine();
				}

				cout << _labelCount[0] << ": " << _labelCount[1] << endl;
				if (_labelCount[1] == 0) {
					iter.NextSlice(); // must be called
					continue;
				}
				cout << "counting done.." << endl;

				time(&now);
				cout << "begin matrix construction " << ctime(&now);

				int count = 0;
				iter.GoToBeginOfSlice();
				int maxRow = 0, maxCol = 0;
				while (!iter.IsAtEndOfSlice()) {
					while (!iter.IsAtEndOfLine()) {
						OrderMapType::IndexType curIdx = iter.GetIndex();
						OrderMapType::PixelType curId = iter.Get();
						LabelMapType::PixelType curLabel = _labelImg->GetPixel(curIdx);
						ImageType::PixelType curVal = _srcImg->GetPixel(curIdx);
						
						NeighborIterator<OrderMapType> nbrIter(_orderImg);	
						nbrIter.GoTo(curIdx);

						double weightSum = 0.0;
						while (!nbrIter.IsAtEnd()) {
							OrderMapType::IndexType nbrIdx = nbrIter.GetIndex();
							OrderMapType::PixelType nbrId = nbrIter.Get();
							LabelMapType::PixelType nbrLabel = _labelImg->GetPixel(nbrIdx);
							ImageType::PixelType nbrVal = _srcImg->GetPixel(nbrIdx);

							double weight = ComputeWeight(curIdx, nbrIdx, curVal, nbrVal);
							if (curId > maxRow) {
								maxRow = curId;
							}
							if (nbrId > maxCol) {
								maxCol = nbrId;
							}
							if (curLabel == 0 && nbrLabel == 0) {
								// update Lu 
								spElement* elm = spGetElement(Lu, curId, nbrId);
								*elm = -weight;
							} else if (curLabel == 0 && nbrLabel > 0) {
								// update B
								spElement* elm = spGetElement(B, curId, nbrId);
								*elm = -weight;
							}
							weightSum += weight;
							++ nbrIter;
						}
						++iter;
					}

					cout << "max Row: " << maxRow << " / max col: " << maxCol << endl;
					iter.NextLine();
				}
				time(&now);
				cout << "building matrix done [" << spElementCount(Lu) << "] " <<  ctime(&now) << endl;

				cout << "Allocating seeds: " << _labelCount[1] << endl;
				spREAL* m = new spREAL[_labelCount[1]];

				cout << "Allocating unseeds: " << _labelCount[0] << endl;
				spREAL* Bm = new spREAL[_labelCount[0]];
				for (int i = 0; i < _labelCount[1]; i++) {
					m[i] = 1;
					Bm[i] = 0;
				}

				spMultiply(B, m, Bm);
				for (int i = 0; i < _labelCount[0]; i++) {
					Bm[i] = -Bm[i];
				}

				spREAL* xs = new spREAL[_labelCount[0]];
				time(&now);
				cout << "Solving matrix ..." << ctime(&now) << endl;
				err = spFactor(Lu);
				if (err >= spFATAL) {
					spErrorMessage(Lu, stderr, "itkRandomWalk");
					return;
				}
				spSolve(Lu, Bm, xs);
				time(&now);
				cout << "Solving matrix ... done" << ctime(&now) << endl;

				iter.GoToBeginOfSlice();
				while (!iter.IsAtEndOfSlice()) {
					while (!iter.IsAtEndOfLine()) {
						OrderMapType::IndexType idx = iter.GetIndex();
						OrderMapType::PixelType id = iter.Get();
						LabelMapType::PixelType label = _labelImg->GetPixel(idx);
						if (label == 0) {
							_probImg->SetPixel(idx, xs[id]);
						}
						++iter;
					}
					iter.NextLine();
				}

				iter.NextSlice();
			}
		}
};

int main(int argc, char* argv[]) {
  if (argc < 4) {
    cout << "usage: " << argv[0] << " input-image input-seed output-prob" << endl;
    return 0;
  }

  int ret = 0;

  ImageType::Pointer srcImg = ReadImageT<ImageType>(argv[1], ret);
  LabelMapType::Pointer seedImg = ReadImageT<LabelMapType>(argv[2], ret);

	RandomWalker rw(srcImg, seedImg);
	rw.Run();

	ImageType::Pointer outImg = rw.GetProbabilityImage(1);	
	WriteImageT<ImageType>(argv[3], outImg);

}
