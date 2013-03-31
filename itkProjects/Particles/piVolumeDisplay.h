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

namespace pi {
    template <class T>
    class VolumeDisplay: public SliceView {
    protected:
        typedef ImageDisplay<T> TDisplay;
        typedef typename T::PixelType TPixel;
        typedef typename T::Pointer TPointer;
        typedef typename T::RegionType TRegion;
        typedef RGBAVolumeType C;
        typedef RGBAVolumeType::Pointer CPointer;
        typedef typename RGBAVolumeType::RegionType CRegion;

        TDisplay *_grayImage;
        CPointer _colorImage;
        std::vector<CPointer> _sliceCache;
        std::vector<void*> _sliceUserData;

        TPixel _windowMin;
        TPixel _windowMax;

        bool _imageModeCache;
        itk::ModifiedTimeType _srcMTime;
        bool _forceUpdate;

    public:
        VolumeDisplay() {
            _grayImage = NULL;
            _srcMTime = 0;
            _imageModeCache = false;
        }

        ~VolumeDisplay() {
        }

        int Count() {
            return _sliceCache.size();
        }

        bool Has(TDisplay* another) {
            return _grayImage == another
                && _grayImage->srcImg == another->srcImg
                && _srcMTime >= another->srcImg->GetMTime();
        }


        void SetDisplay(TDisplay* display) {
            _grayImage = display;
            _srcMTime = display->srcImg->GetMTime();
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

        bool UpdateSlice(SliceDirectionEnum dir, bool navigationImage = false) {
            if (_grayImage == NULL) {
                return false;
            }

            TPixel newMin = _grayImage->histogram.rangeMin;
            TPixel newMax = _grayImage->histogram.rangeMax;

            if (!_forceUpdate && dir == _sliceDirection && _windowMin == newMin && _windowMax == newMax && _imageModeCache == navigationImage) {
                return false;
            }

            _sliceDirection = dir;
            _windowMin = newMin;
            _windowMax = newMax;
            _imageModeCache = navigationImage;
            _forceUpdate = false;

            _sliceCache.clear();
            _sliceUserData.clear();

            if (navigationImage) {
                _colorImage = VolumeDisplay<T>::ConvertToColor(_grayImage->navigationImg, _grayImage->histogram.rangeMin, _grayImage->histogram.rangeMax);
            } else {
                _colorImage = VolumeDisplay<T>::ConvertToColor(_grayImage->srcImg, _grayImage->histogram.rangeMin, _grayImage->histogram.rangeMax);
            }

            CRegion region = _colorImage->GetBufferedRegion();
            typedef itk::ExtractImageFilter<C, C> ExtractFilter;
            const int nSlices = region.GetSize(_sliceDirection);
            for (int i = 0; i < nSlices; i++) {
                CPointer slice = SliceDisplay<RGBAVolumeType>::ExtractSlice(_colorImage, i, _sliceDirection);
                if (i == 0) {
                    SliceView view = SliceDisplay<RGBAVolumeType>::GetSliceView(slice);
                    _width = view.Width();
                    _height = view.Height();
                }
                _sliceCache.push_back(slice);
            }

            _sliceUserData.resize(_sliceCache.size());
            for (int i = 0; i < _sliceUserData.size(); i++) {
                _sliceUserData[i] = NULL;
            }
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