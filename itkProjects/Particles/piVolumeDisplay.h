//
//  piVolumeView.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/31/13.
//
//

#ifndef ParticleGuidedRegistration_piVolumeView_h
#define ParticleGuidedRegistration_piVolumeView_h

#include "piImageSlice.h"
#include <itkExtractImageFilter.h>
#include "itkScalarToARGBColormapImageFilter.h"
#include <algorithm>
#include <QMutex>
#include <QMutexLocker>

static QMutex __mutex;

namespace pi {
    template <class T>
    class VolumeDisplay: public SliceView {
    protected:
        typedef ImageDisplay<T> TDisplay;
        typedef typename ImageDisplay<T>::Pointer TDisplayPointer;
        typedef typename T::PixelType TPixel;
        typedef typename T::Pointer TPointer;
        typedef typename T::RegionType TRegion;
        typedef RGBAVolumeType C;
        typedef RGBAVolumeType::Pointer CPointer;
        typedef typename RGBAVolumeType::RegionType CRegion;

        TDisplayPointer _grayImage;
        CPointer _colorImage;
        std::vector<CPointer> _sliceCache;
        std::vector<void*> _sliceUserData;

        TPixel _windowMin;
        TPixel _windowMax;

        bool _imageModeCache;
        itk::ModifiedTimeType _srcMTime;
        bool _forceUpdate;

        SliceDirectionEnum _processingDirection;
        TPixel _processingMin;
        TPixel _processingMax;
        bool _underProcessing;

    public:
        VolumeDisplay() {
            _grayImage = NULL;
            _srcMTime = 0;
            _imageModeCache = false;
            _underProcessing = false;
        }

        ~VolumeDisplay() {
        }

        int Count() {
            return _sliceCache.size();
        }

        bool Has(TDisplayPointer another) {
            // FIXME: test pointer error
            return _grayImage == another && _srcMTime >= _grayImage->GetMTime();
        }

        void SetDisplay(TDisplayPointer display) {
            _grayImage = display;
            _srcMTime = display->GetMTime();
            _forceUpdate = true;
        }

        void SetSliceData(int i, void* pointer) {
            _sliceUserData[i] = pointer;
        }

        template <class S>
        S* GetSliceData(int i) {
            if (i < 0 || i >= _sliceCache.size()) {
                return (S*) NULL;
            } else {
                return (S*) _sliceUserData[i];
            }
        }

        unsigned char* GetColorImageBuffer(int i) {
            if (i < 0 || i >= _sliceCache.size()) {
                return NULL;
            }
            return (uchar*) _sliceCache[i]->GetBufferPointer();
        }

        void operator()() {
            _processingDirection = _sliceDirection;
            CPointer colorImage = VolumeDisplay<T>::ConvertToColor(_grayImage->GetSourceImage(), _grayImage->GetHistogram().rangeMin, _grayImage->GetHistogram().rangeMax);

            CRegion region = colorImage->GetBufferedRegion();
            const int nSlices = region.GetSize(_sliceDirection);

            std::vector<CPointer> sliceCache;
            typedef itk::ExtractImageFilter<C, C> ExtractFilter;
            for (int i = 0; i < nSlices; i++) {
                CPointer slice = SliceDisplay<RGBAVolumeType>::ExtractSlice(colorImage, i, _processingDirection);
                if (i == 0) {
                    SliceView view = SliceDisplay<RGBAVolumeType>::GetSliceView(slice);
                    _width = view.Width();
                    _height = view.Height();
                }
                sliceCache.push_back(slice);
            }

            QMutexLocker lock(&__mutex);
            _colorImage = colorImage;
            _sliceCache = sliceCache;
            

        }

        bool UpdateSlice(SliceDirectionEnum dir, bool navigationImage = false) {
            if (_grayImage.IsNull()) {
                return false;
            }

            TPixel newMin = _grayImage->GetHistogram().rangeMin;
            TPixel newMax = _grayImage->GetHistogram().rangeMax;

            bool changeSliceDirection = (dir != _sliceDirection);
            bool changeIntensityWindow = (_windowMin != newMin || _windowMax != newMax);

            if (!_forceUpdate && (!changeSliceDirection) && (!changeIntensityWindow) && _imageModeCache == navigationImage) {
                return false;
            }

            _sliceDirection = dir;
            _windowMin = newMin;
            _windowMax = newMax;
            _imageModeCache = navigationImage;
            _forceUpdate = false;

            CRegion region = _grayImage->GetSourceImage()->GetBufferedRegion();
            const int nSlices = region.GetSize(_sliceDirection);

            if (changeSliceDirection || _sliceUserData.size() == 0) {
                _sliceUserData.clear();
                _sliceUserData.resize(nSlices);
                for (int i = 0; i < nSlices; i++) {
                    _sliceUserData[i] = NULL;
                }
            }

            this->operator()();
            return true;
        }

        static CPointer ConvertToColor(TPointer grayImage, TPixel viewMin, TPixel viewMax, itk::ColormapEnumType color = itk::Grey) {
            // convert to color image
            // assume that the composite image ranges from 0 to 65535 (in short range)
            typedef itk::ScalarToARGBColormapImageFilter<T, RGBAVolumeType> ColorFilterType;
            typename ColorFilterType::Pointer colorFilter = ColorFilterType::New();
            colorFilter->SetInput(grayImage);
            colorFilter->UseManualScalingOn();
            colorFilter->SetMinimumValue(viewMin);
            colorFilter->SetMaximumValue(viewMax);
            colorFilter->Update();
            CPointer rgbImage = colorFilter->GetOutput();
            rgbImage->DisconnectPipeline();
            return rgbImage;
        }
        
        TPointer& operator[](int i) {
            return _sliceCache[i];
        }

    };
}

#endif