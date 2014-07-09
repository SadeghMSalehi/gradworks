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

#include "itkObject.h"
#include "itkSmartPointer.h"
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
    
    typedef pi::RealImage AIRImage;
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

    template <class T>
    class ImageDisplayCollection;
    
    template<class T>
    class ImageDisplay: public itk::Object {
    public:
        typedef ImageDisplay Self;
        typedef Object Superclass;
        typedef ImageDisplayCollection<T> CollectionType;
        typedef SliceDisplay<T> SliceType;
        typedef itk::SmartPointer<Self> Pointer;
        typedef itk::SmartPointer<const Self> ConstPointer;

        itkTypeMacro(ImageDisplay, Superclass);
        itkNewMacro(ImageDisplay);

        typedef itk::AffineTransform<double> TransformType;
        typedef TransformType::OutputVectorType VectorType;
        typedef SliceDisplay<T> GridType;

        typedef typename T::Pointer ImagePointer;
        typedef typename T::PointType ImagePoint;
        typedef typename T::SpacingType ImageSpacing;
        typedef ImageHistogram<T> HistogramType;

        itkGetMacro(SourceImage, ImagePointer);
        itkGetMacro(Origin, ImagePoint);
        itkGetMacro(Spacing, ImageSpacing);
        itkGetMacro(FileName, std::string);
        itkGetMacro(AffineTransform, TransformType::Pointer);

        itkSetMacro(FileName, std::string);
        itkSetMacro(AffineTransform, TransformType::Pointer);

        itkGetMacro(ParentCollection, CollectionType*);

    protected:
        ImageDisplay() {
            m_ParentCollection = NULL;
        }
        virtual ~ImageDisplay() {
        }

    private:
        ImageDisplay(const Self& other);
        void operator=(const Self& other);

        CollectionType* m_ParentCollection;
        ImagePointer m_SourceImage;
        ImagePoint m_Origin;
        ImageSpacing m_Spacing;

        std::string m_FileName;
        itk::RescaledImage<T> navigationImg;
        TransformType::Pointer m_AffineTransform;

        typename T::PointType _currentOrigin;
        VectorType _currentTranslation;

        // (grid x transform)
        std::vector<SliceType> _sliceDisplayCache;

        HistogramType histogram;
        ImageDisplayProperty displayProperty;

    public:
        HistogramType& GetHistogram() {
            return histogram;
        }

        void SetParentCollection(CollectionType* parent) {
            m_ParentCollection = parent;
            InvalidateCaches();
        }

        void SetAffineTransform(vtkMatrix4x4* mat) {
            if (m_SourceImage.IsNull()) {
                return;
            }
            TransformType::ParametersType params = m_AffineTransform->GetParameters();
            for (int i = 0; i < DIMENSIONS; i++) {
                for (int j = 0; j < DIMENSIONS; j++) {
                    params[3 * i + j] = mat->GetElement(i, j);
                }
            }
            m_AffineTransform->SetParameters(params);
            InvalidateCaches();
        }

        void SetAffineTransform(itk::TransformBase::Pointer baseTransform) {
            if (m_SourceImage.IsNull()) {
                return;
            }
            m_AffineTransform->SetParameters(baseTransform->GetParameters());
            m_AffineTransform->SetFixedParameters(baseTransform->GetFixedParameters());
            InvalidateCaches();
        }


        void SetAffineParameters(TransformType::ParametersType params) {
            if (params.GetSize() == 12) {
                m_AffineTransform->SetParameters(params);
            }
            InvalidateCaches();
        }

        VectorType GetAffineTranslation() {
            return m_AffineTransform->GetTranslation();
        }
        
        void SetAffineTranslation() {
            _currentTranslation = m_AffineTransform->GetTranslation();
        }

        void SetAffineTranslation(int idx, float x, float y) {
            GridType& refGrid = m_ParentCollection->GetDisplay(idx);
            TransformType::OutputVectorType translation = _currentTranslation;
            if (refGrid.Direction() == IJ) {
                translation[0] += x * m_Spacing[0];
                translation[1] += y * m_Spacing[1];
            } else if (refGrid.Direction() == JK) {
                translation[1] += x * m_Spacing[1];
                translation[2] += y * m_Spacing[2];
            } else if (refGrid.Direction() == KI) {
                translation[0] += x * m_Spacing[0];
                translation[2] += y * m_Spacing[2];
            }
            m_AffineTransform->SetTranslation(translation);
            InvalidateCaches();
        }

        void SetAffineTranslation(VectorType tx) {
            m_AffineTransform->SetTranslation(tx);
            InvalidateCaches();
        }

        void SetAffineCenter(typename T::PointType center) {
            m_AffineTransform->SetCenter(center);
            InvalidateCaches();
        }

        // 2d slice related functions
        void SetSourceImage(typename T::Pointer img) {
            if (img.IsNull()) {
                cout << "Empty source image" << endl;
                return;
            }

            m_SourceImage = img;
            m_Origin = img->GetOrigin();
            m_Spacing = img->GetSpacing();
            
            histogram.SetImage(m_SourceImage);

            TransformType::ParametersType params;
            params.SetSize(12);
            params.Fill(0);
            for (int i = 0; i < DIMENSIONS; i++) {
                params[3 * i + i] = 1;
            }
            
            m_AffineTransform = TransformType::New();
            m_AffineTransform->SetParameters(params);

            typename T::PointType center;
            typename T::SizeType imgSz = m_SourceImage->GetBufferedRegion().GetSize();
            RealIndex centerIdx;
            fordim (k) {
                centerIdx[k] = imgSz[k]/2.0;
            }
            m_SourceImage->TransformContinuousIndexToPhysicalPoint(centerIdx, center);
            m_AffineTransform->SetCenter(center);

            InvalidateCaches();
        }

        void SetWindowRange(typename T::PixelType m1, typename T::PixelType m2) {
            displayProperty.windowMin = m1;
            displayProperty.windowMax = m2;
        }

        void SetColorMap(itk::ColormapEnumType map) {
            displayProperty.colorMap = map;
        }

        typename T::Pointer GetDisplayImage(int j = 0) {
            if (_sliceDisplayCache.size() != m_ParentCollection->GetDisplayCount()) {
                _sliceDisplayCache.resize(m_ParentCollection->GetDisplayCount());
            }
            if (_sliceDisplayCache[j].IsNull()) {
                Resample(j);
            }
            return _sliceDisplayCache[j];
        }

        void InvalidateCaches(int i = -1) {
            if (i == -1) {
                if (m_ParentCollection != NULL) {
                    _sliceDisplayCache = std::vector<SliceType>(m_ParentCollection->GetDisplayCount());
                }
            } else {
                _sliceDisplayCache[i] = SliceType();
            }
        }

        typename T::Pointer Resample3D(typename T::Pointer resampleGrid, int interpolator = 0) {
            typedef itk::ResampleImageFilter<AIRImage, AIRImage> ResampleFilter;
            typename ResampleFilter::Pointer resampleFilter = ResampleFilter::New();
            resampleFilter->SetInput(m_SourceImage);
            resampleFilter->SetReferenceImage(resampleGrid);
            resampleFilter->UseReferenceImageOn();
            if (interpolator == 1) {
                typedef itk::BSplineInterpolateImageFunction<T> BspFuncType;
                resampleFilter->SetInterpolator(BspFuncType::New());
            }
            if (m_AffineTransform.IsNotNull()) {
                typename ResampleFilter::TransformType* transformInput = dynamic_cast<typename ResampleFilter::TransformType*>(m_AffineTransform.GetPointer());
                resampleFilter->SetTransform(transformInput);
            }
            resampleFilter->Update();
            return resampleFilter->GetOutput();
        }

    private:
        void Resample(int j) {
            if (m_ParentCollection == NULL || m_ParentCollection->GetDisplay(j).IsNull()) {
                cout << "No sampling grid" << endl;
                return;
            }
            if (m_SourceImage.IsNull()) {
                cout << "No source image" << endl;
                return;
            }
            if (_sliceDisplayCache.size() < m_ParentCollection->GetDisplayCount()) {
                _sliceDisplayCache.resize(m_ParentCollection->GetDisplayCount());
            }
            if (_sliceDisplayCache[j].IsNotNull()) {
                return;
            }
            _sliceDisplayCache[j] = Resample3D(m_ParentCollection->GetDisplay(j));
        }
    };


    template<class T>
    class ImageDisplayCollection {
    public:
        typedef SliceDisplay<T> SliceType;
        typedef ImageDisplay<T> ImageDisplayType;
        typedef typename ImageDisplayType::Pointer ImageDisplayPointer;

    private:
        std::vector<ImageDisplayPointer> _imageDisplays;
        std::vector<SliceType> _displayArray;
        typename T::PointType _referenceCenter;
        typename T::RegionType _referenceRegion;
        typedef typename T::PixelType TPixel;

        void SetCenterOfRotation() {
            if (_imageDisplays.size() > 0) {
                typename T::IndexType idx1 = _referenceRegion.GetIndex();
                typename T::IndexType idx2 = _referenceRegion.GetUpperIndex();
                RealIndex centerIdx;
                for (int i = 0; i < 3; i++) {
                    centerIdx[i] = (idx1[i] + idx2[i]) / 2.0;
                }
                _imageDisplays[0]->GetSourceImage()->TransformContinuousIndexToPhysicalPoint(centerIdx, _referenceCenter);
                for (int i = 0; i < _imageDisplays.size(); i++) {
                    _imageDisplays[i].SetAffineCenter(_referenceCenter);
                    _imageDisplays[i].InvalidateCaches();
                }
            }
        }

    public:
        ImageDisplayCollection() {
            _imageDisplays.reserve(2);
        }

        bool IsEmpty() {
            return _imageDisplays.size() == 0;
        }

        void Reset() {
            _imageDisplays.clear();
            _displayArray.clear();
            _referenceCenter.Fill(0);
            _referenceRegion = typename  T::RegionType();
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
        void SetWindowRange(int n, TPixel min, TPixel max) {
            _imageDisplays[n].SetWindowRange(min, max);
        }

        // image management
        ImageDisplayPointer AddImage(typename T::Pointer srcImg, std::string fileName = "") {
            ImageDisplayPointer newImage = ImageDisplayType::New();
            newImage->SetSourceImage(srcImg);
            newImage->SetFileName(fileName);
            newImage->SetParentCollection(this);
            if (_imageDisplays.size() == 0) {
                _referenceRegion = srcImg->GetBufferedRegion();
            }
            if (_imageDisplays.size() > 1) {
                newImage->SetAffineCenter(_referenceCenter);
            }
            _imageDisplays.push_back(newImage);
            return newImage;
        }


        bool SetImage(int n, typename T::Pointer img) {
            if (!IsValidId(n)) {
                return false;
            }
            _imageDisplays[n]->SetSourceImage(img);
            return true;
        }

        ImageDisplayPointer& at(int i) {
            return _imageDisplays[i];
        }

        ImageDisplayPointer& operator[](int i) {
            return _imageDisplays[i];
        }

        ImageDisplayPointer& GetReference() {
            return _imageDisplays[0];
        }

        int GetReferenceSize(SliceDirectionEnum dir) {
            return _referenceRegion.GetSize(dir);
        }
        
        ImageDisplayPointer& GetLast() {
            return _imageDisplays.back();
        }

        // recreate a slice grid
        void SetSliceDisplay(SliceDirectionEnum axis, int index) {
            if (_imageDisplays.size() == 0) {
                return;
            }
            typename T::Pointer srcImg = _imageDisplays[0]->GetSourceImage();
            if (srcImg.IsNull()) {
                cout << "Emtpy reference grid" << endl;
                return;
            }

            typename T::Pointer sliceImg = SliceDisplay<T>::ExtractSlice(srcImg, index, axis);
            if (sliceImg.IsNull()) {
                return;
            }
            SetDisplay(sliceImg);
        }

        void SetDisplay(typename T::Pointer grid) {
            _displayArray.clear();
            _displayArray.push_back(SliceType(grid));
            for (int i = 0; i < _imageDisplays.size(); i++) {
                _imageDisplays[i]->InvalidateCaches();
            }
        }

        SliceType& GetDisplay(int i) {
            return _displayArray[i];
        }

        int GetDisplayCount() {
            return _displayArray.size();
        }
    };

    typedef std::vector<RGBAVolumeType::Pointer> RGBAImageVector;
    typedef std::vector<AIRImage::Pointer> AIRImageVector;
    typedef SliceDisplay<AIRLabel> AIRLabelSlice;
    typedef SliceDisplay<AIRImage> AIRImageSlice;
    typedef ImageDisplay<AIRImage>::Pointer AIRImageDisplay;
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
