//
//  piImageSlice.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 2/23/13.
//
//

#ifndef __ParticleGuidedRegistration__piImageSlice__
#define __ParticleGuidedRegistration__piImageSlice__

#include <iostream>
#include <algorithm>

#include "piImageDef.h"
#include "piImageHistogram.h"

#include "itkRGBAPixel.h"
#include "itkScalarToARGBColormapImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkTransform.h"
#include "itkAffineTransform.h"
#include "vtkMatrix4x4.h"
#include "itkBSplineInterpolateImageFunction.h"

#include "QPixmap"

class QGraphicsPixmapItem;

namespace itk {
    template <class T>
    class PointerCache {
    private:
        bool _cached;

    public:
        typename T::Pointer Cache;

        PointerCache() {
            _cached = false;
        }

        PointerCache(typename T::Pointer data) {
            Cache = data;
            _cached = true;
        }

        ~PointerCache() {
            Cache = NULL;
        }

        bool HasSamePointer(typename T::Pointer another) {
            return Cache.GetPointer() == another.GetPointer();
        }

        bool IsValid() {
            return Cache.IsNotNull() && _cached;
        }

        void Invalidate() {
            _cached = false;
            Cache = NULL;
        }

        PointerCache& operator=(typename T::Pointer data) {
            if (data.IsNotNull()) {
                Cache = data;
                _cached = true;
            } else {
                Cache = NULL;
                _cached = false;
            }
            return (*this);
        }
    };
}

/**
 * this is a class representing a slice in QT gui
 */
namespace pi {
    typedef itk::RGBAPixel<unsigned char> RGBAPixel;
    typedef itk::Image<RGBAPixel, 2> RGBAImageType;
    typedef itk::Image<RGBAPixel, 3> RGBAVolumeType;
    
    typedef pi::LabelImage AIRImage;
    typedef AIRImage::PixelType AIRPixel;

    class ImageDisplayProperty {
    public:
        AIRPixel windowMin;
        AIRPixel windowMax;
        itk::ColormapEnumType colorMap;

        ImageDisplayProperty() {
            colorMap = itk::Grey;
        }
    };

    enum SliceDirectionEnum { IJ = 2, JK = 0, KI = 1, Unknown = -1 };

    template<class T>
    class ImageResampleGrid {
    private:
        typename T::Pointer GridImg;
        int _width;
        int _height;
        SliceDirectionEnum SliceDirection;

    public:
        ImageResampleGrid(typename T::Pointer gridImg) {
            GridImg = gridImg;
            
            typename T::RegionType region = gridImg->GetBufferedRegion();
            if (region.GetSize(0) == 1) {
                SliceDirection = JK;
                _width = region.GetSize(1);
                _height = region.GetSize(2);
            } else if (region.GetSize(1) == 1) {
                SliceDirection = KI;
                _width = region.GetSize(0);
                _height = region.GetSize(2);
            } else if (region.GetSize(2) == 1) {
                _width = region.GetSize(0);
                _height = region.GetSize(1);
                SliceDirection = IJ;
            }
        }

        SliceDirectionEnum Direction() {
            return SliceDirection;
        }
        
        typename T::ConstPointer GetConstImage() {
            return typename T::ConstPointer(GridImg);
        }
        
        operator typename T::Pointer() {
            return GridImg;
        }

        operator const T*() {
            return GridImg.GetPointer();
        }

        bool IsNull() {
            return GridImg.IsNull();
        }

        bool IsNotNull() {
            return GridImg.IsNotNull();
        }

        inline int Width() {
            return _width;
        }

        inline int Height() {
            return _height;
        }
    };

    template<class T>
    class ImageDisplay {
    public:
        typedef itk::AffineTransform<double> TransformType;
        typedef TransformType::OutputVectorType VectorType;
        typedef ImageResampleGrid<T> GridType;

    private:
        TransformType::Pointer _affineTransform;
        typename T::PointType _currentOrigin;
        VectorType _currentTranslation;

        // (grid x transform)
        typedef itk::PointerCache<T> ImageCache;
        std::vector<ImageCache> _resampledImageCache;

    public:
        typename T::Pointer srcImg;
        typename T::PointType srcOrigin;
        typename T::SpacingType srcSpacing;

        // how to determine bin-size automatically?
        ImageHistogram<T> histogram;
        ImageDisplayProperty displayProperty;

        // set from outside
        std::vector<GridType>* _resampleGridVector;

    public:
        ImageDisplay(typename T::Pointer img) {
            SetImage(img);
            _resampleGridVector = NULL;
        }

        ~ImageDisplay() {
            _resampleGridVector = NULL;
            _resampledImageCache.clear();
        }

        void SetGridVector(std::vector<GridType>* grids) {
            _resampleGridVector = grids;
            _resampledImageCache.clear();
            _resampledImageCache.resize(grids->size());
        }

        void SetAffineTransform(vtkMatrix4x4* mat) {
            if (srcImg.IsNull()) {
                return;
            }
            TransformType::ParametersType params = _affineTransform->GetParameters();
            for (int i = 0; i < DIMENSIONS; i++) {
                for (int j = 0; j < DIMENSIONS; j++) {
                    params[3 * i + j] = mat->GetElement(i, j);
                }
            }
            _affineTransform->SetParameters(params);
            InvalidateCaches();
        }

        void SetAffineTransform(TransformType::TransformBase::Pointer baseTransform) {
            if (srcImg.IsNull()) {
                return;
            }
            _affineTransform->SetParameters(baseTransform->GetParameters());
            _affineTransform->SetFixedParameters(baseTransform->GetFixedParameters());
            InvalidateCaches();
        }

        TransformType::Pointer GetAffineTransform() {
            return _affineTransform;
        }

        void SetAffineParameters(TransformType::ParametersType params) {
            if (params.GetSize() == 12) {
                _affineTransform->SetParameters(params);
            }
            InvalidateCaches();
        }

        VectorType GetAffineTranslation() {
            return _affineTransform->GetTranslation();
        }
        
        void SetAffineTranslation() {
            _currentTranslation = _affineTransform->GetTranslation();
        }

        void SetAffineTranslation(int idx, float x, float y) {
            GridType& refGrid = _resampleGridVector->at(idx);
            TransformType::OutputVectorType translation = _currentTranslation;
            if (refGrid.Direction() == IJ) {
                translation[0] += x * srcSpacing[0];
                translation[1] += y * srcSpacing[1];
            } else if (refGrid.Direction() == JK) {
                translation[1] += x * srcSpacing[1];
                translation[2] += y * srcSpacing[2];
            } else if (refGrid.Direction() == KI) {
                translation[0] += x * srcSpacing[0];
                translation[2] += y * srcSpacing[2];
            }
            _affineTransform->SetTranslation(translation);
            InvalidateCaches();
        }

        void SetAffineTranslation(VectorType tx) {
            _affineTransform->SetTranslation(tx);
            InvalidateCaches();
        }

        void SetAffineCenter(typename T::PointType center) {
            _affineTransform->SetCenter(center);
            InvalidateCaches();
        }

        // 2d slice related functions
        void SetImage(typename T::Pointer img) {
            if (img.IsNull()) {
                cout << "Empty source image" << endl;
                return;
            }

            srcImg = img;
            srcOrigin = img->GetOrigin();
            srcSpacing = img->GetSpacing();
            
            histogram.SetImage(srcImg);
//            histogram.Compute();

            TransformType::ParametersType params;
            params.SetSize(12);
            params.Fill(0);
            for (int i = 0; i < DIMENSIONS; i++) {
                params[3 * i + i] = 1;
            }
            
            _affineTransform = TransformType::New();
            _affineTransform->SetParameters(params);

            typename T::PointType center;
            typename T::SizeType imgSz = srcImg->GetBufferedRegion().GetSize();
            RealIndex centerIdx;
            fordim (k) {
                centerIdx[k] = imgSz[k]/2.0;
            }
            srcImg->TransformContinuousIndexToPhysicalPoint(centerIdx, center);
            _affineTransform->SetCenter(center);

            InvalidateCaches();
        }

        void SetWindowRange(typename T::PixelType m1, typename T::PixelType m2) {
            displayProperty.windowMin = m1;
            displayProperty.windowMax = m2;
        }

        void SetColorMap(itk::ColormapEnumType map) {
            displayProperty.colorMap = map;
        }

        SliceDirectionEnum GetResampleDirection(int j) {
            if (_resampleGridVector == NULL || j < 0 || _resampleGridVector->size() <= j || _resampleGridVector->at(j).IsNull()) {
                return Unknown;
            }
            return _resampleGridVector->at(j).Direction();
        }

        typename T::Pointer GetResampled(int j = 0) {
            if (j < 0) {
                cout << "Wrong sampling index" << endl;
                return typename T::Pointer();
            }
            if (_resampleGridVector == NULL || _resampleGridVector->size() <= j || _resampleGridVector->at(j).IsNull()) {
                cout << "no sampling grid" << endl;
                return typename T::Pointer();
            }
            if (_resampleGridVector->size() != _resampledImageCache.size()) {
                _resampledImageCache.resize(_resampleGridVector->size());
            }
            if (!_resampledImageCache[j].IsValid()) {
                Resample(j);
            }
            return _resampledImageCache[j].Cache;
        }

        void InvalidateCaches(int i = -1) {
            if (i == -1) {
                for (int j = 0; j < _resampledImageCache.size(); j++) {
                    _resampledImageCache[j].Invalidate();
                }
            } else {
                _resampledImageCache[i].Invalidate();
            }
        }

        typename T::Pointer Resample3D(typename T::Pointer resampleGrid, int interpolator = 0) {
            typedef itk::ResampleImageFilter<AIRImage, AIRImage> ResampleFilter;
            typename ResampleFilter::Pointer resampleFilter = ResampleFilter::New();
            resampleFilter->SetInput(srcImg);
            resampleFilter->SetReferenceImage(resampleGrid);
            resampleFilter->UseReferenceImageOn();
            if (interpolator == 1) {
                typedef itk::BSplineInterpolateImageFunction<T> BspFuncType;
                resampleFilter->SetInterpolator(BspFuncType::New());
            }
            if (_affineTransform.IsNotNull()) {
                typename ResampleFilter::TransformType* transformInput = dynamic_cast<typename ResampleFilter::TransformType*>(_affineTransform.GetPointer());
                resampleFilter->SetTransform(transformInput);
            }
            resampleFilter->Update();
            return resampleFilter->GetOutput();
        }

    private:
        void Resample(int j) {
            if (_resampleGridVector->at(j).IsNull()) {
                cout << "No sampling grid" << endl;
                return;
            }
            if (srcImg.IsNull()) {
                cout << "No source image" << endl;
                return;
            }
            if (_resampledImageCache[j].IsValid()) {
                return;
            }
            _resampledImageCache[j] = Resample3D(_resampleGridVector->at(j));
        }
    };


    template<class T>
    class ImageDisplayCollection {
    public:
        typedef ImageResampleGrid<T> GridType;
        typedef ImageDisplay<T> ImageDisplayType;
        typedef std::vector<ImageDisplayType> ImageDisplayVector;

    private:
        int _referenceId;
        typename T::PointType _referenceCenter;
        ImageDisplayVector _imageDisplays;
        std::vector<GridType> _resampleGrids;

        void SetCenterOfRotation() {
            if (IsValidId(_referenceId)) {
                typename T::Pointer refImg = _imageDisplays[_referenceId].srcImg;
                typename T::RegionType refRegion = refImg->GetBufferedRegion();
                typename T::IndexType idx1 = refRegion.GetIndex();
                typename T::IndexType idx2 = refRegion.GetUpperIndex();
                RealIndex centerIdx;
                for (int i = 0; i < 3; i++) {
                    centerIdx[i] = (idx1[i] + idx2[i]) / 2.0;
                }
                refImg->TransformContinuousIndexToPhysicalPoint(centerIdx, _referenceCenter);

                for (int i = 0; i < _imageDisplays.size(); i++) {
                    _imageDisplays[i].SetAffineCenter(_referenceCenter);
                    _imageDisplays[i].InvalidateCaches();
                }
            }
        }

    public:
        ImageDisplayCollection(): _referenceId(-1) {
            _imageDisplays.reserve(2);
        }

        void Reset() {
            _imageDisplays.clear();
            _resampleGrids.clear();
            _referenceId = -1;
            _referenceCenter.Fill(0);
        }

        bool IsValidId(int n) {
            if (n >= 0 && n < _imageDisplays.size()) {
                return true;
            }
            return false;
        }

        int Count() {
            return _imageDisplays.size();
        }
        
        // may be changed to a piece-wise function
        void SetWindowRange(int n, AIRPixel min, AIRPixel max) {
            _imageDisplays[n].SetWindowRange(min, max);
        }

        // image management
        ImageDisplayType& AddImage(typename T::Pointer srcImg) {
            ImageDisplayType newImage(srcImg);
            _imageDisplays.push_back(ImageDisplayType(srcImg));

            ImageDisplayType& lastImage = _imageDisplays.back();
            lastImage.SetGridVector(&_resampleGrids);

            if (IsValidId(_referenceId)) {
                lastImage.SetAffineCenter(_referenceCenter);
            }
            return lastImage;
        }
//
//        void DeleteImage(int n) {
//            if (n == _referenceId && _imageDisplays.size() > 1) {
//                SetReferenceId(0);
//            }
//            if (n >= 0 && n < _imageDisplays.size()) {
//                _imageDisplays.erase(_imageDisplays.begin() + n);
//            }
//            if (_imageDisplays.size() == 0) {
//                _referenceId = -1;
//            }
//        }

        bool SetImage(int n, typename T::Pointer img) {
            if (IsValidId(n)) {
                _imageDisplays[n].SetImage(img);
            } else if (n == Count()) {
                _imageDisplays.push_back(ImageDisplayType(img));
            } else {
                return false;
            }
            _imageDisplays[n].SetGridVector(&_resampleGrids);
            if (n == _referenceId) {
                SetReferenceId(n);
            }
            return true;
        }

        ImageDisplayType* at(int i) {
            return &(_imageDisplays[i]);
        }
        GridType& GetGrid(int i) {
            return _resampleGrids[i];
        }
        ImageDisplayType& operator[](int i) {
            return _imageDisplays[i];
        }
        ImageDisplayType& GetReference() {
            return _imageDisplays[_referenceId];
        }
        typename T::Pointer GetReferenceGrid() {
            if (IsValidId(_referenceId)) {
                return _imageDisplays[_referenceId].srcImg;
            } else {
                return typename T::Pointer();
            }
        }
        int GetReferenceSize(SliceDirectionEnum dir) {
            return GetReferenceGrid()->GetBufferedRegion().GetSize(dir);
        }
        ImageDisplayType& GetLast() {
            return _imageDisplays.back();
        }

        // what should happen if the reference is changed?
        void SetReferenceId(int n) {
            if (n >= 0 && n < _imageDisplays.size()) {
                if (_referenceId != n) {
                    _referenceId = n;
                    SetCenterOfRotation();
                }
            }
        }

        // recreate a slice grid
        void SetSliceGrid(SliceDirectionEnum axis, int index) {
            typename T::Pointer srcImg = GetReferenceGrid();
            if (srcImg.IsNull()) {
                cout << "Emtpy reference grid" << endl;
                return;
            }

            typename T::RegionType resampleRegion = srcImg->GetBufferedRegion();
            typename T::IndexType idx1 = resampleRegion.GetIndex();
            typename T::IndexType idx2 = resampleRegion.GetUpperIndex();

            idx1[axis] = index;
            idx2[axis] = index;
            resampleRegion.SetIndex(idx1);
            resampleRegion.SetUpperIndex(idx2);

            typedef itk::ExtractImageFilter<T, T> ExtractFilterType;
            typename ExtractFilterType::Pointer filter = ExtractFilterType::New();
            filter->SetInput(srcImg);
            filter->SetExtractionRegion(resampleRegion);
            filter->Update();
            
            SetResampleGrid(filter->GetOutput());
        }

        void SetResampleGrid(typename T::Pointer grid) {
            _resampleGrids.clear();
            _resampleGrids.push_back(GridType(grid));
            for (int i = 0; i < _imageDisplays.size(); i++) {
                _imageDisplays[i].InvalidateCaches();
            }
        }
    };

    typedef std::vector<RGBAVolumeType::Pointer> RGBAImageVector;
    typedef std::vector<AIRImage::Pointer> AIRImageVector;
    typedef ImageDisplay<AIRImage> AIRDisplayImage;
    typedef std::vector<AIRDisplayImage> AIRDISplayVector;
    typedef ImageDisplayCollection<AIRImage> AIRDisplayCollection;

    
    template<class T>
    class ImagePixmap {
    public:
        int alpha;
        QGraphicsPixmapItem* pixmapCache;
        
        ImagePixmap() :
        alpha(255), pixmapCache(NULL) {
        }
        
        void SetImage(typename T::Pointer label) {
            labelImage = label;
            typedef itk::ScalarToARGBColormapImageFilter<T, RGBAImageType> ScalarToRGBFilter;
            typename ScalarToRGBFilter::Pointer rgbFilter = ScalarToRGBFilter::New();
            rgbFilter->SetInput(labelImage);
            rgbFilter->UseManualScalingOff();
            rgbFilter->SetAlphaValue(alpha);
            rgbFilter->Update();
            rgbaImage = rgbFilter->GetOutput();
        }
        
        typename T::Pointer GetImage() {
            return labelImage;
        }
        
        QPixmap GetPixmap() {
            if (rgbaImage.IsNull()) {
                return QPixmap();
            }
            RGBAImageType::SizeType bitmapSz = rgbaImage->GetBufferedRegion().GetSize();
            QImage qImg = QImage((unsigned char*) rgbaImage->GetBufferPointer(), bitmapSz[0], bitmapSz[1], QImage::Format_ARGB32);
            return QPixmap::fromImage(qImg);
        }
    private:
        typename T::Pointer labelImage;
        RGBAImageType::Pointer rgbaImage;
    };
}
#endif /* defined(__ParticleGuidedRegistration__piImageSlice__) */
