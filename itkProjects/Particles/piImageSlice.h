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

    typedef itk::ScalarToARGBColormapImageFilter<RealImage, RGBAVolumeType> ColorFilterType;

    class ImageDisplayProperty {
    public:
        DataReal windowMin;
        DataReal windowMax;
        ColorFilterType::ColormapEnumType colorMap;

        ImageDisplayProperty() {
            colorMap = ColorFilterType::Grey;
        }
    };

    enum SliceDirectionEnum { IJ = 2, JK = 0, KI = 1 };

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
        typedef ImageResampleGrid<T> GridType;
        
        typename T::Pointer srcImg;
        typename T::PointType srcOrigin;
        typename T::SpacingType srcSpacing;

        // how to determine bin-size automatically?
        ImageHistogram<T> histogram;
        ImageDisplayProperty displayProperty;
        
        typedef itk::AffineTransform<double> TransformType;
        TransformType::Pointer affineTransform;

        // added by user
        std::vector<GridType>* _resampleGridVector;
    private:
        typename T::PointType _currentOrigin;

        // (grid x transform)
        typedef itk::PointerCache<T> ImageCache;
        std::vector<ImageCache> _resampledImageCache;


    public:
        ImageDisplay(typename T::Pointer img) {
            SetImage(img);
            _resampleGridVector = NULL;
        }

        ~ImageDisplay() {
            _resampleGridVector = NULL;
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
            TransformType::ParametersType params = affineTransform->GetParameters();
            for (int i = 0; i < DIMENSIONS; i++) {
                for (int j = 0; j < DIMENSIONS; j++) {
                    params[3 * i + j] = mat->GetElement(i, j);
                }
            }
            affineTransform->SetParameters(params);
            InvalidateCaches();
        }

        void SetCenterOfRotation(typename T::PointType center) {
            TransformType::ParametersType centerParams = affineTransform->GetFixedParameters();
            for (int i = 0; i < 3; i++) {
                centerParams[i] = center[i];
            }
            affineTransform->SetFixedParameters(centerParams);
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
            
            affineTransform = TransformType::New();
            affineTransform->SetParameters(params);
            InvalidateCaches();
        }


        /*
        void MoveOrigin(double x, double y, double z) {
            typename T::PointType origin = srcImg->GetOrigin();
            origin[0] += x * srcSpacing[0];
            origin[1] += y * srcSpacing[1];
            origin[2] += z * srcSpacing[2];
            srcImg->SetOrigin(origin);
            InvalidateCaches();
        }
         */

        void SetOrigin() {
            _currentOrigin = srcImg->GetOrigin();
        }

        void SetOrigin(double x, double y, double z) {
            typename T::PointType origin = _currentOrigin;
            origin[0] += x * srcSpacing[0];
            origin[1] += y * srcSpacing[1];
            origin[2] += z * srcSpacing[2];
            srcImg->SetOrigin(origin);
            InvalidateCaches();
        }

        void SetOriginFromResampling(int idx, float x, float y) {
            GridType& refGrid = _resampleGridVector->at(idx);
            typename T::PointType origin = _currentOrigin;
            if (refGrid.Direction() == IJ) {
                origin[0] += x * srcSpacing[0];
                origin[1] += y * srcSpacing[1];
            } else if (refGrid.Direction() == JK) {
                origin[1] += x * srcSpacing[1];
                origin[2] += y * srcSpacing[2];
            } else if (refGrid.Direction() == KI) {
                origin[0] += x * srcSpacing[0];
                origin[2] += y * srcSpacing[2];
            }
            srcImg->SetOrigin(origin);
            InvalidateCaches();
        }

        void SetOrigin(int axis, double t) {
            typename T::PointType origin = _currentOrigin;
            origin[axis] = t * srcSpacing[axis];
            srcImg->SetOrigin(origin);
            InvalidateCaches();
        }

        void SetWindowRange(typename T::PixelType m1, typename T::PixelType m2) {
            displayProperty.windowMin = m1;
            displayProperty.windowMax = m2;
        }

        void SetColorMap(ColorFilterType::ColormapEnumType map) {
            displayProperty.colorMap = map;
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
            typedef itk::ResampleImageFilter<RealImage, RealImage> ResampleFilter;
            typename ResampleFilter::Pointer resampleFilter = ResampleFilter::New();
            resampleFilter->SetInput(srcImg);
            resampleFilter->SetReferenceImage(_resampleGridVector->at(j));
            resampleFilter->UseReferenceImageOn();
            if (affineTransform.IsNotNull()) {
                typename ResampleFilter::TransformType* transformInput = dynamic_cast<typename ResampleFilter::TransformType*>(affineTransform.GetPointer());
                resampleFilter->SetTransform(transformInput);
            }
            resampleFilter->Update();
            _resampledImageCache[j] = resampleFilter->GetOutput();
//            cout << __FILE__ << ":" << __LINE__ << " Resampling done" << endl;
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
                    _imageDisplays[i].SetCenterOfRotation(_referenceCenter);
                    _imageDisplays[i].InvalidateCaches();
                }
            }
        }

    public:
        ImageDisplayCollection(): _referenceId(-1) {
            _imageDisplays.reserve(10);
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
        void SetWindowRange(int n, pi::DataReal min, pi::DataReal max) {
            _imageDisplays[n].SetWindowRange(min, max);
        }

        // image management
        ImageDisplayType& AddImage(pi::RealImage::Pointer srcImg) {
            ImageDisplayType newImage(srcImg);
            _imageDisplays.push_back(ImageDisplayType(srcImg));

            ImageDisplayType& lastImage = _imageDisplays.back();
            lastImage.SetGridVector(&_resampleGrids);

            if (IsValidId(_referenceId)) {
                lastImage.SetCenterOfRotation(_referenceCenter);
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
//        void SetImage(int n, pi::RealImage::Pointer img) {
//            if (n >= 0 && n < _imageDisplays.size()) {
//                _imageDisplays[n].SetImage(img);
//                if (n == _referenceId) {
//                    SetCenterOfRotation();
//                }
//            }
//        }
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
            return _imageDisplays[_referenceId].srcImg;
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
            if (_referenceId < 0 || _imageDisplays[_referenceId].srcImg.IsNull()) {
                return;
            }

            typename T::Pointer srcImg = _imageDisplays[_referenceId].srcImg;
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
