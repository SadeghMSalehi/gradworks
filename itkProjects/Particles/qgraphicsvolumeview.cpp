//
//  qgraphicsvolumeview.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/27/13.
//
//

#include "qgraphicsvolumeview.h"

#include <itkResampleImageFilter.h>
#include <itkExtractImageFilter.h>
#include <QPixmap>
#include <QImage>
#include <QGraphicsPixmapItem>

using namespace pi;

QGraphicsVolumeView::QGraphicsVolumeView(QWidget* parent): QGraphicsView(parent) {
    setScene(&_scene);

    _thumbsWidth = 128;
    _columnCount = 4;

    _displayId = 0;
    _displayReference = true;
    _manualIntensityScaling = false;
    _directionCache = Unknown;

    _airImages = NULL;
}

void QGraphicsVolumeView::setDisplayCollection(pi::AIRDisplayCollection *images) {
    if (images == NULL) {
        return;
    }

    if (_airImages == images) {
        return;
    }

    _airImages = images;
    _scene.clear();
    _sliceCache.clear();
    _displayImages.clear();
    _directionCache = Unknown;
    _volumeCache = NULL;
    
    updateDisplay();
}

bool QGraphicsVolumeView::checkSliceCache() {
    AIRDisplayImage& refImg = _airImages->GetReference();
    SliceDirectionEnum dir = refImg.GetResampleDirection(0);

    // fire only when the slice direction is changed
    if (_directionCache == dir) {
        return true;
    }

    _directionCache = dir;
    return false;
}

bool QGraphicsVolumeView::updateSource() {
    if (_directionCache == Unknown) {
        return false;
    }
    AIRImage::RegionType region = _volumeCache->GetBufferedRegion();

    _sliceCache.clear();
    _sliceCache.reserve(region.GetSize(_directionCache));

    typedef itk::ExtractImageFilter<AIRImage, AIRImage> ExtractFilter;

    const int nSlices = region.GetSize(_directionCache);
    for (int i = 0; i < nSlices; i++) {
        ExtractFilter::Pointer extract = ExtractFilter::New();
        extract->SetInput(_volumeCache);
        region.SetIndex(_directionCache, i);
        region.SetSize(_directionCache, 1);
        _volumeCache->GetBufferedRegion().Print(cout);
        region.Print(cout);
        extract->SetExtractionRegion(region);
        extract->Update();
        _sliceCache.push_back(extract->GetOutput());
        _sliceCache.back()->DisconnectPipeline();
    }
    return true;
}

void QGraphicsVolumeView::updateDisplay() {
    if (this->parentWidget()->isHidden()) {
        return;
    }
    
    if (_airImages == NULL) {
        return;
    }

    if (!checkVolumeCache()) {
        return;
    }

    if (!checkSliceCache()) {
        if (!updateSource()) {
            return;
        }
    }

    // delete all graphics items
    _displayImages.clear();
    _scene.clear();

    // fire when intensity changed
    AIRImageVector::iterator sliceIter;

    int colPos = 0;
    int rowPos = 0;

    AIRImage::SizeType volumeSize = _volumeCache->GetBufferedRegion().GetSize();
    int w = volumeSize[0];
    int h = volumeSize[1];
    if (_directionCache == JK) {
        w = volumeSize[1];
        h = volumeSize[2];
    } else if (_directionCache == KI) {
        w = volumeSize[0];
        h = volumeSize[2];
    }

    AIRDisplayImage& refImg = _airImages->GetReference();
    for (sliceIter = _sliceCache.begin(); sliceIter != _sliceCache.end(); sliceIter++) {
        typedef itk::ScalarToARGBColormapImageFilter<AIRImage, RGBAVolumeType> ColorFilterType;
        ColorFilterType::Pointer colorFilter = ColorFilterType::New();
        colorFilter->SetInput((*sliceIter));
        colorFilter->UseManualScalingOn();
        colorFilter->SetMinimumValue(refImg.histogram.rangeMin);
        colorFilter->SetMaximumValue(refImg.histogram.rangeMax);
        colorFilter->Update();
        RGBAVolumeType::Pointer rgbImg = colorFilter->GetOutput();
        rgbImg->DisconnectPipeline();
        _displayImages.push_back(rgbImg);
        QGraphicsPixmapItem* item = _scene.addPixmap(QPixmap::fromImage(QImage((unsigned char*)rgbImg->GetBufferPointer(), w, h, QImage::Format_ARGB32)));
        item->translate(w*colPos, h*rowPos);
        if (++colPos % _columnCount == 0) {
            colPos = 0;
            rowPos ++;
        }
    }
}

bool QGraphicsVolumeView::checkVolumeCache() {
    AIRDisplayImage& refImg = _airImages->GetReference();
    if (refImg.srcImg.IsNull()) {
        _volumeCache = NULL;
        return false;
    }
    if (_volumeCache != refImg.srcImg) {
        _volumeCache = NULL;
    }

    AIRImage::SpacingType spacing = refImg.srcSpacing;
    AIRImage::SizeType size = refImg.srcImg->GetBufferedRegion().GetSize();

    double resampleFactor = size[0] / _thumbsWidth;
    fordim(k) {
        spacing[k] = spacing[k] * resampleFactor;
        size[k] = size[k] / resampleFactor;
    }

    typedef itk::ResampleImageFilter<AIRImage, AIRImage> ResampleFilter;
    ResampleFilter::Pointer resampler = ResampleFilter::New();
    resampler->SetInput(refImg.srcImg);
    resampler->SetOutputParametersFromImage(refImg.srcImg.GetPointer());
    resampler->SetOutputSpacing(spacing);
    resampler->SetSize(size);
    resampler->Update();
    _volumeCache = resampler->GetOutput();
    _volumeCache->DisconnectPipeline();

    return _volumeCache.IsNotNull();
}