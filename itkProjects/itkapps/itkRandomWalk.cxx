#include "itkImageCommon.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageSliceIteratorWithIndex.h"
#include "itkConstNeighborhoodIterator.h"
#include "math.h"
#include "stdio.h"
#include "time.h"
#include "spMatrix.h"
#include "iostream"
#include "vector"


using namespace std;

typedef itk::Image<short,3> ImageType;
typedef itk::Image<double,3> DoubleImageType;
typedef itk::ImageRegionIteratorWithIndex<ImageType> IteratorType;
typedef itk::ImageRegionIteratorWithIndex<DoubleImageType> DoubleIteratorType;
typedef itk::ImageSliceIteratorWithIndex<ImageType> SliceIteratorType;
typedef itk::ConstNeighborhoodIterator<ImageType> NeighborhoodIterator;
typedef itk::ImageSliceIteratorWithIndex<DoubleImageType> DoubleImageSliceIteratorType;

void test() {
	spMatrix A;
	spError err;

	A = spCreate(10, 0, &err);
	spClear(A);

	for (int i = 1; i <= 10; i++) {
		spElement* e = spGetElement(A, i, i);
		*e = i;
	}

	for (int i = 1; i <= 10; i++) {
		for (int j = 1; j <= 10; j++) {
			spElement* e = NULL;
			if (NULL != (e = spFindElement(A, i, j))) {
				cout << *e << " ";
			} else {
				cout << 0 << " ";
			}
		}
		cout << endl;
	}

	cout << "number of elements: ";
	cout << spElementCount(A) << endl;
}

double ComputeWeight(DoubleImageType::IndexType curIdx, DoubleImageType::IndexType nbrIdx, DoubleImageType::PixelType curVal, DoubleImageType::PixelType nbrVal) {
	double weight = exp(-0.5*(curVal - nbrVal)*(curVal - nbrVal));
	return weight;
}

spREAL ProcessNeighbor(ImageType::Pointer seedImg, ImageType::Pointer srcImg, 
												ImageType::IndexType curIdx, ImageType::IndexType nbrIdx, 
												int curSub, int nbrSub, spMatrix& Lu, spMatrix& B) {

	ImageType::PixelType curSeed = seedImg->GetPixel(curIdx);
	ImageType::PixelType nbrSeed = seedImg->GetPixel(nbrIdx);

	ImageType::PixelType curVal = srcImg->GetPixel(curIdx);
	ImageType::PixelType nbrVal = srcImg->GetPixel(nbrIdx);

	spElement* elm = NULL;
	if (curSeed == 1 && nbrSeed == 1) {
		return 0;
	}

	if (curSeed == 0 && nbrSeed == 0) {
		elm = spGetElement(Lu, curSub, nbrSub);
	} else if (curSeed > 0 && nbrSeed == 0) {
		elm = spGetElement(B, nbrSub, curSub);
	} else if (curSeed == 0 && nbrSeed > 0) {
		elm = spGetElement(B, curSub, nbrSub);
	}

	double weight = ComputeWeight(curIdx, nbrIdx, curVal, nbrVal);
	*elm = -weight;

	return *elm;
}


spREAL* BuildMatrix(ImageType::Pointer srcImg, ImageType::Pointer seedImg) {
	time_t now;
	ImageType::RegionType seedRegion = seedImg->GetBufferedRegion();
	ImageType::SizeType imgSize = seedRegion.GetSize();

	int curSub = 0;
	vector<spREAL> seeds(100000);
	vector<spREAL> unseeds(100000);

	spMatrix Lu, B;
	spError err;
	Lu = spCreate(1000000, 0, &err);
	B = spCreate(1000000, 0, &err);

	int seedsCount = 0;
	int unseedsCount = 0;

	time(&now);
	cout << "Processing neighborhoods ..." << ctime(&now) << endl;

	NeighborhoodIterator::RadiusType radius;
	for (int i = 0; i < 3; i++) { 
		radius[i] = 1;
	}

	NeighborhoodIterator iter(radius, seedImg, seedRegion);
	for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
		ImageType::PixelType curLabel = iter.GetCenterPixel();

		if (curLabel > 0) {
			seedsCount ++;
			seeds.push_back(curLabel);
		} else {
			unseedsCount++;
		}

		ImageType::IndexType curIdx = iter.GetIndex();
		double weightSum = 0.0;
		for (int i = 0; i < 9; i++) {
			ImageType::IndexType nbrIdx = iter.GetIndex(i);
			if (seedRegion.IsInside(nbrIdx)) {
				int nbrSub = nbrIdx[0] + nbrIdx[1] * imgSize[0] + nbrIdx[2] * imgSize[0] * imgSize[1];
				double weight = ProcessNeighbor(seedImg, srcImg, curIdx, nbrIdx, curSub, nbrSub, Lu, B);
				weightSum = weightSum + abs(weight);
			} 
		}

		spElement* elm = spGetElement(Lu, curSub, curSub);
		*elm = weightSum;

		++curSub;
	}

	time(&now);
	cout << "Processing neighborhoods ... done : " << ctime(&now) << endl; 

	spREAL* m = &seeds[0];
	spREAL* Bm = new spREAL[seedsCount];

	cout << "RHS multiplication " << endl;
	spMultiply(B, m, Bm);

	for (int i = 0; i < seedsCount; i++) {
		Bm[i] = -Bm[i];
	}

	spREAL* xs = new spREAL[unseedsCount];
	time(&now);
	cout << "Solving matrix ..." << ctime(&now) << endl; 
	spSolve(Lu, Bm, xs);
	time(&now);
	cout << "Solving matrix ... done" << ctime(&now) << endl; 

	delete[] Bm;
	return xs;
}

DoubleImageType::Pointer BuildProbabilityImage(ImageType::Pointer seedImg, spREAL* xs) {
	ImageType::SizeType szImg = seedImg->GetBufferedRegion().GetSize();
	DoubleImageType::Pointer probImg = NewImageT<DoubleImageType>(szImg[0], szImg[1], szImg[2]);

	IteratorType iter(seedImg, seedImg->GetBufferedRegion());
	DoubleIteratorType outIter(probImg, probImg->GetBufferedRegion());

	int xsIdx = 0;
	for (; !iter.IsAtEnd(); ++iter, ++outIter) {
		if (iter.Get() > 0) {
			outIter.Set(0);
		} else {
			outIter.Set(xs[xsIdx]);
			xsIdx++;
		}
	}

	return probImg;
}

void IterateNeighbors(ImageType::IndexType idx) {

}


template <class TImage>
class ImageAccessor {
	public:
		ImageAccessor(typename TImage::Pointer img) {
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
			_xyOffsets = { {-1,0,-1}, {0,0,-1}, {1,0,-1}, {-1,0,0}, {1,0,0}, {-1,0,1}, {0,0,1}, {1,0,1} };
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

void UpdateMatrix() {
	int sub = accr.SliceIndex2Sub(idx);
	ImageType::IndexType sub2idx = accr.Sub2SliceIndex(sliceIdx, sub);

	NeighborIterator<ImageType> nbrIter(seedImg);
	nbrIter.GoTo(idx);

	double weightSum = 0.0;
	while (!nbrIter.IsAtEnd()) {
		ImageType::IndexType nbrIdx = nbrIter.GetIndex();
		ImageType::PixelType nbrVal = nbrIter.Get();
		int nbrsub = accr.SliceIndex2Sub(nbrIdx);

		double weight = ComputeWeight(idx, nbrIdx, val, nbrVal);
		if (val == 1 && nbrVal == 0) {
			spElement* elm = spGetElement(B, nbrsub, sub);
			*elm = -weight;
		} else if (val == 0 && nbrVal == 1) {
			spElement* elm = spGetElement(B, sub, nbrsub);
			*elm = -weight;
		} else if (val == 0 && nbrVal == 0) {
			spElement* elm = spGetElement(Lu, sub, nbrsub);
			*elm = -weight;
		}

		weightSum = weightSum + weight;
		++nbrIter;
	}

	spElement* elm = spGetElement(Lu, sub, sub);
	*elm = weightSum;
}

int main(int argc, char* argv[]) {
	if (argc < 4) {
		cout << "usage: " << argv[0] << " input-image input-seed output-mask output-prob" << endl;
		return 0;
	}

	int ret = 0;
	time_t now;

	DoubleImageType::Pointer srcImg = ReadImageT<DoubleImageType>(argv[1], ret);
	ImageType::Pointer seedImg = ReadImageT<ImageType>(argv[2], ret);

	ImageType::SizeType srcSz = srcImg->GetBufferedRegion().GetSize();
	DoubleImageType::Pointer probImg = NewImageT<DoubleImageType>(srcSz[0], srcSz[1], srcSz[2]);
	DoubleImageSliceIteratorType outIter(probImg, probImg->GetBufferedRegion());

	SliceIteratorType iter(seedImg, seedImg->GetBufferedRegion());
	iter.SetFirstDirection(0);
	iter.SetSecondDirection(2);

	ImageAccessor<ImageType> accr(seedImg);

	spMatrix Lu, B;
	spError err;
	Lu = spCreate(1000000, 0, &err);
	B  = spCreate(1000000, 0, &err);

	int sliceIdx = 0;
	while (!iter.IsAtEnd()) {
		int seedsCount = 0;
		int unseedsCount = 0;

		std::vector<ImageType::IndexType> label0Index;
		std::vector<ImageType::IndexType> label1Index;


		while (!iter.IsAtEndOfSlice()) {
			while (!iter.IsAtEndOfLine()) {
				ImageType::PixelType val = iter.Get();
				ImageType::IndexType idx = iter.GetIndex();

				if (val == 0) {
					label0Index.push_back(idx);
				} else if (val == 1) {
					label1Index.push_back(idx);
				}
				iter.NextLine();
			}

			if (label1Index.size() > 0) {
				UpdateMatrix();

			}

			iter.NextSlice();
		}


				++iter;
			}
			iter.NextLine();
		}

		cout << "Allocating seeds: " << seedsCount << endl;
		spREAL* m = new spREAL[seedsCount];

		cout << "Allocating unseeds: " << unseedsCount << endl;
		spREAL* Bm = new spREAL[unseedsCount];
		for (int i = 0; i < seedsCount; i++) {
			m[i] = 1;
			Bm[i] = 0;
		}
		
		spMultiply(B, m, Bm);
		for (int i = 0; i < unseedsCount; i++) {
			Bm[i] = -Bm[i];
		}

		spREAL* xs = new spREAL[unseedsCount];
		time(&now);
		cout << "Solving matrix ..." << ctime(&now) << endl; 
		spSolve(Lu, Bm, xs);
		time(&now);
		cout << "Solving matrix ... done" << ctime(&now) << endl; 


		for (int i = 0; i < (int) unseedsIndex.size(); i++) {
			ImageType::IndexType unseedsPos = unseedsIndex[i];
			probImg->SetPixel(unseedsPos, xs[i]);
		}

		outIter.NextSlice();
		iter.NextSlice();

		cout << "Processing slice: " << ++sliceIdx << endl;
	}
	/*
	ImageType::SizeType srcSize = srcImg->GetBufferedRegion().GetSize();
	cout << srcSize << endl;

	for (int i = 0; i < srcSize[1]; i++) {
		spREAL* xs = BuildMatrix2D(srcImg, seedImg, i);
		DoubleImageType::Pointer probImg = BuildProbabilityImage2D(seedImg, xs);
	}
	*/

	WriteImageT<DoubleImageType>(argv[2], probImg);
}
