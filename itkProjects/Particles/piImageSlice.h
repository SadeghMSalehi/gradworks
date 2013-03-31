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
#include "piImageIO.h"

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
#pragma mark -
#pragma Templated Class Definitions

#pragma mark -
#pragma Templated Class Definitions

    template <class T>
    class RescaledImage {
    private:
        typename T::Pointer _img;
        float _rescaleFactor;
    public:
        RescaledImage() {
            _rescaleFactor = 1;
        }

        bool IsNull() {
            return _img.IsNull();
        }

        bool IsNotNull() {
            return _img.IsNotNull();
        }
        
        ulong GetMTime() {
            return _img->GetMTime();
        }

        operator typename T::Pointer() {
            return _img;
        }

        typename T::SizeType GetSize() {
            return _img->GetBufferedRegion().GetSize();
        }

        int GetSize(int dir) {
            return _img->GetBufferedRegion().GetSize(dir);
        }

        int GetOriginalIndex(int idx) {
            return idx * _rescaleFactor;
        }

        int GetIndexFromOriginal(int idx) {
            return idx / _rescaleFactor;
        }

        void SetImage(typename T::Pointer img, float rescale = 1) {
            if (img.IsNull() || (_img.IsNotNull() && _img->GetMTime() > img->GetMTime())) {
                return;
            }
            _rescaleFactor = rescale;

            typename T::SpacingType spacing = img->GetSpacing();
            typename T::SizeType size = img->GetBufferedRegion().GetSize();

            fordim(k) {
                spacing[k] = spacing[k] * _rescaleFactor;
                size[k] = size[k] / _rescaleFactor;
            }

            typedef itk::ResampleImageFilter<T,T> ResampleFilter;
            typename ResampleFilter::Pointer resampler = ResampleFilter::New();
            resampler->SetInput(img);
            resampler->SetOutputParametersFromImage(img.GetPointer());
            resampler->SetOutputSpacing(spacing);
            resampler->SetSize(size);
            resampler->Update();

            _img = resampler->GetOutput();
            _img->DisconnectPipeline();
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
    typedef itk::Image<uchar,3> AIRLabel;
    typedef AIRImage::PixelType AIRPixel;
    typedef AIRLabel::PixelType AIRClass;

    class ImageDisplayProperty {
    public:
        AIRPixel windowMin;
        AIRPixel windowMax;
        itk::ColormapEnumType colorMap;

        ImageDisplayProperty() {
            colorMap = itk::Grey;
        }
    };

    class SliceView {
    protected:
        int _width;
        int _height;
        int _index;
        SliceDirectionEnum _sliceDirection;
    public:
        SliceView(): _width(0), _height(0), _index(-1), _sliceDirection(Unknown) {}
        SliceView(int w, int h, int i, SliceDirectionEnum d): _width(w), _height(h), _index(i), _sliceDirection(d) {}
        virtual ~SliceView() {}
        inline int Width() { return _width; }
        inline int Height() { return _height; }
        inline int Index() { return _index; }
        inline SliceDirectionEnum Direction() { return _sliceDirection; }
    };

    template<class T>
    class SliceDisplay: public SliceView {
    protected:
        typedef typename T::PixelType TPixel;
        typedef typename T::Pointer TPointer;
        typedef typename T::RegionType TRegion;
        TPointer _srcImg;
        int _depth;

    public:
        SliceDisplay() : _depth(0) {
        }
        
        SliceDisplay(typename T::Pointer gridImg): SliceView(0,0,-1,Unknown), _depth(0) {
            SetImage(gridImg);
        }

        virtual ~SliceDisplay() {
        }

        void SetImage(TPointer img) {
            _srcImg = img;
            if (_srcImg.IsNull()) {
                cout << "source image is null " << __FILE__ << ":" << __LINE__ << endl;
                return;
            }
            SliceView view = SliceDisplay<T>::GetSliceView(_srcImg);
            _width = view.Width();
            _height = view.Height();
            _index = view.Index();
            _sliceDirection = view.Direction();
            if (_sliceDirection == Unknown) {
                _depth = _srcImg->GetBufferedRegion().GetSize(2);
            }
        }

        int Index() {
            return _index;
        }

        typename T::ConstPointer GetConstImage() {
            return typename T::ConstPointer(_srcImg);
        }

        typename T::Pointer GetImage() {
            return _srcImg;
        }
        
        operator typename T::Pointer() {
            return _srcImg;
        }

        operator const T*() {
            return _srcImg.GetPointer();
        }

        bool IsNull() {
            return _srcImg.IsNull();
        }

        bool IsNotNull() {
            return _srcImg.IsNotNull();
        }
        inline int Depth() {
            return _depth;
        }

        typename T::PixelType* GetBufferPointer() {
            if (_srcImg.IsNull()) {
                return NULL;
            }
            return _srcImg->GetBufferPointer();
        }

        typename T::RegionType GetRegion() {
            if (_srcImg.IsNull()) {
                return typename T::RegionType();
            }
            return _srcImg->GetBufferedRegion();
        }
        
        typename T::RegionType GetRegion(int left, int top, int width, int height) {
            if (_srcImg.IsNull()) {
                return typename T::RegionType();
            }
            typename T::RegionType region =  _srcImg->GetBufferedRegion();
            switch (_sliceDirection) {
                case IJ:
                    region.SetIndex(0, left);
                    region.SetIndex(1, top);
                    region.SetSize(0, width);
                    region.SetSize(1, height);
                    break;
                case JK:
                    region.SetIndex(1, left);
                    region.SetIndex(2, top);
                    region.SetSize(1, width);
                    region.SetSize(2, height);
                    break;
                case KI:
                    region.SetIndex(0, left);
                    region.SetIndex(2, top);
                    region.SetSize(0, width);
                    region.SetSize(2, height);
                    break;
                default:
                    break;
            }
            return region;
        }

        void FillBuffer(typename T::PixelType p) {
            if (_srcImg.IsNull()) {
                return;
            }
            _srcImg->FillBuffer(p);
        }


        static SliceView GetSliceView(typename T::Pointer img) {
            typename T::RegionType region = img->GetBufferedRegion();
            if (region.GetSize(0) == 1) {
                return SliceView(region.GetSize(1), region.GetSize(2), region.GetIndex(0), JK);
            } else if (region.GetSize(1) == 1) {
                return SliceView(region.GetSize(0), region.GetSize(2), region.GetIndex(1), KI);
            } else if (region.GetSize(2) == 1) {
                return SliceView(region.GetSize(0), region.GetSize(1), region.GetIndex(2), IJ);
            } else {
                return SliceView(region.GetSize(0), region.GetSize(2), region.GetIndex(0), Unknown);
            }
        }

        static SliceDisplay<T> ExtractSlice(typename T::Pointer srcImg, int idx, SliceDirectionEnum dir) {
            typename T::Pointer emptyImg;
            if (srcImg.IsNull() || dir == Unknown) {
                cout << "source is null or direction unknown:" << __FILE__ << ":" << __LINE__ << endl;
                return emptyImg;
            }

            typename T::RegionType sliceRegion = srcImg->GetBufferedRegion();
            typename T::SizeType srcSize = sliceRegion.GetSize();
            if (srcSize[dir] <= idx || idx < 0) {
                cout << "slice index is wrong: " << idx << " >= " << srcSize[dir] << __FILE__ << ":" << __LINE__ << endl;
                return emptyImg;
            }
            sliceRegion.SetIndex(dir, idx);
            sliceRegion.SetSize(dir,1);

            typedef itk::ExtractImageFilter<T, T> ExtractFilterType;
            typename ExtractFilterType::Pointer filter = ExtractFilterType::New();
            filter->SetInput(srcImg);
            filter->SetExtractionRegion(sliceRegion);
            filter->Update();
            typename T::Pointer sliceImg = filter->GetOutput();
            sliceImg->DisconnectPipeline();
            return SliceDisplay<T>(sliceImg);
        }
    };

    template<class T>
    class ImageDisplay {
    public:
        typedef itk::AffineTransform<double> TransformType;
        typedef TransformType::OutputVectorType VectorType;
        typedef SliceDisplay<T> GridType;

    private:
        TransformType::Pointer _affineTransform;
        typename T::PointType _currentOrigin;
        VectorType _currentTranslation;

        // (grid x transform)
        std::vector<typename T::Pointer> _resampledImageCache;

    public:
        typename T::Pointer srcImg;
        typename T::PointType srcOrigin;
        typename T::SpacingType srcSpacing;
        std::string fileName;

        itk::RescaledImage<T> navigationImg;

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

        itk::RescaledImage<T>& GetNavigationImage() {
            if (navigationImg.IsNull() && srcImg->GetMTime() > navigationImg.GetMTime()) {
                navigationImg.SetImage(srcImg, 2);
            }
            return navigationImg;
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
            if (!_resampledImageCache[j].IsNull()) {
                Resample(j);
            }
            return _resampledImageCache[j];
        }

        void InvalidateCaches(int i = -1) {
            if (i == -1) {
                for (int j = 0; j < _resampledImageCache.size(); j++) {
                    _resampledImageCache[j] = NULL;
                }
            } else {
                _resampledImageCache[i] = NULL;
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
            if (_resampledImageCache[j].IsNull()) {
                return;
            }
            _resampledImageCache[j] = Resample3D(_resampleGridVector->at(j));
        }
    };


    


    template<class T>
    class ImageDisplayCollection {
    public:
        typedef SliceDisplay<T> GridType;
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

        int GetReferenceId() {
            return _referenceId;
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
        ImageDisplayType& AddImage(typename T::Pointer srcImg, std::string fileName = "") {
            ImageDisplayType newImage(srcImg);
            newImage.fileName = fileName;
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
        void SetReferenceSlice(SliceDirectionEnum axis, int index) {
            typename T::Pointer srcImg = GetReferenceGrid();
            if (srcImg.IsNull()) {
                cout << "Emtpy reference grid" << endl;
                return;
            }

            typename T::Pointer sliceImg = SliceDisplay<T>::ExtractSlice(srcImg, index, axis);
            if (sliceImg.IsNull()) {
                return;
            }
            SetReferenceSlice(sliceImg);
        }

        void SetReferenceSlice(typename T::Pointer grid) {
            _resampleGrids.clear();
            _resampleGrids.push_back(GridType(grid));
            for (int i = 0; i < _imageDisplays.size(); i++) {
                _imageDisplays[i].InvalidateCaches();
            }
        }
    };

    typedef std::vector<RGBAVolumeType::Pointer> RGBAImageVector;
    typedef std::vector<AIRImage::Pointer> AIRImageVector;
    typedef SliceDisplay<AIRLabel> AIRLabelSlice;
    typedef SliceDisplay<AIRImage> AIRImageSlice;
    typedef ImageDisplay<AIRImage> AIRImageDisplay;
    typedef std::vector<AIRImageDisplay> AIRDisplayVector;
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

    extern ImageIO<AIRImage> __airImageIO;
    extern ImageIO<AIRLabel> __airLabelIO;
}
#endif /* defined(__ParticleGuidedRegistration__piImageSlice__) */
