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
#include <QMouseEvent>
#include <QList>
#include <QKeyEvent>
#include <QGraphicsRectItem>

using namespace pi;

QGraphicsVolumeView::QGraphicsVolumeView(QWidget* parent): QGraphicsView(parent) {
    setScene(&_scene);

    _thumbsWidth = 128;
    _columnCount = 10000;

    _displayId = 0;
    _displayReference = true;
    _manualIntensityScaling = false;
    _directionCache = Unknown;

    _airImages = NULL;

    setDragMode(QGraphicsView::RubberBandDrag);
//    setBackgroundBrush(QBrush(QPixmap::fromImage(QImage(QString::fromUtf8(":/Icons/Images/backgroundPattern.jpg")))));
}

void QGraphicsVolumeView::setDisplayCollection(pi::AIRDisplayCollection *images) {
    if (images == NULL) {
        return;
    }

    if (_airImages == images) {
        return;
    }
    clear();

    _airImages = images;

    updateDisplay();
}

void QGraphicsVolumeView::clear() {
    _scene.clear();
    _sliceCache.clear();
    _slicePixmaps.clear();
    _displayImages.clear();
    _workingSet.clear();
    _directionCache = Unknown;
    _volumeCache = NULL;
}

std::vector<int> QGraphicsVolumeView::getWorkingSet() {
    std::vector<int> workingSet;
    QIntSet::ConstIterator iter = _workingSet.constBegin();
    for (; iter != _workingSet.constEnd(); iter++) {
        int sliceIdx = *iter * _rescaleFactor;
        int nextSliceIdx = sliceIdx + 1 * _rescaleFactor;
        for (int i = sliceIdx; i < nextSliceIdx; i++) {
            workingSet.push_back(i);
        }
    }
    return workingSet;
}

bool QGraphicsVolumeView::checkVolumeCache() {
    AIRDisplayImage& refImg = _airImages->GetReference();
    if (refImg.srcImg.IsNull()) {
        _volumeCache = NULL;
        return false;
    }
    if (_volumeCache != refImg.srcImg) {
        clear();
    }

    AIRImage::SpacingType spacing = refImg.srcSpacing;
    AIRImage::SizeType size = refImg.srcImg->GetBufferedRegion().GetSize();

    _rescaleFactor = 1;
    if (_thumbsWidth != 0) {
        _rescaleFactor = size[0] / _thumbsWidth;
    }

    if (std::abs(_rescaleFactor - 1) < 0.1) {
        _volumeCache = refImg.srcImg;
        return _volumeCache.IsNotNull();
    }

    fordim(k) {
        spacing[k] = spacing[k] * _rescaleFactor;
        size[k] = size[k] / _rescaleFactor;
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
        extract->SetExtractionRegion(region);
        extract->Update();
        _sliceCache.push_back(extract->GetOutput());
        _sliceCache.back()->DisconnectPipeline();
    }
    return true;
}

void QGraphicsVolumeView::updateDisplay() {
    if (this->isHidden()) {
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
    _slicePixmaps.clear();

    // fire when intensity changed
    AIRImageVector::iterator sliceIter;

    int colPos = 0;

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

    _slicePixmaps.clear();

    AIRDisplayImage& refImg = _airImages->GetReference();
    int sliceIdx = 0;
    for (sliceIter = _sliceCache.begin(); sliceIter != _sliceCache.end(); sliceIter++, sliceIdx++) {
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
        item->setFlags(QGraphicsItem::ItemIsSelectable);
        item->translate(colPos, 0);
        item->setData(SliceIndex, QVariant(sliceIdx));
        _slicePixmaps.push_back(item);

        if (_workingSet.contains(sliceIdx)) {
            highlightSlice(item);
        }

        QGraphicsTextItem* text = _scene.addText(QString("%1").arg(sliceIdx), QFont("Courier", 24));
        text->setPos(colPos+3, 3);
        text->setZValue(1);
        text->setDefaultTextColor(Qt::yellow);


        colPos += (w + 1);
    }
}

void QGraphicsVolumeView::highlightSlice(QGraphicsPixmapItem *sliceItem) {
    QGraphicsRectItem* rectItem = new QGraphicsRectItem(sliceItem->boundingRect(), sliceItem);
    rectItem->setZValue(1);
    rectItem->setPen(QPen(Qt::yellow, 3));
}


void QGraphicsVolumeView::keyReleaseEvent(QKeyEvent* key) {
    if (key->text() == "=") {
        scale(1.1, 1.1);
        this->setMaximumHeight(this->maximumHeight()*1.2);
    } else if (key->text() == "-") {
        scale(1/1.1,1/1.1);
        this->setMaximumHeight(this->maximumHeight()/1.1);
    } else if (key->text() == "0") {
        this->resetTransform();
    }
}

void QGraphicsVolumeView::createWorkingSet() {
    QList<QGraphicsItem*> selections = _scene.selectedItems();
    QList<QGraphicsItem*>::ConstIterator itemIter = selections.constBegin();
    for (;itemIter != selections.constEnd(); itemIter++) {
        QGraphicsPixmapItem* sliceItem = dynamic_cast<QGraphicsPixmapItem*>(*itemIter);
        if (sliceItem == NULL) {
            continue;
        }
        int sliceIdx = sliceItem->data(SliceIndex).value<int>();
        if (_workingSet.contains(sliceIdx)) {
            _scene.removeItem(sliceItem->childItems()[0]);
            _workingSet.remove(sliceIdx);
        } else {
            _workingSet.insert(sliceIdx);
            highlightSlice(sliceItem);
        }
    }
    this->setInteractive(true);
}

void QGraphicsVolumeView::clearWorkingSet() {
    QIntSet::ConstIterator iter = _workingSet.constBegin();
    for (; iter != _workingSet.constEnd(); iter++) {
        _scene.removeItem(_slicePixmaps[*iter]->childItems()[0]);
    }
    _workingSet.clear();
    this->setInteractive(true);
}

void QGraphicsVolumeView::mousePressEvent(QMouseEvent* event) {
    if (event->buttons() & Qt::RightButton) {
        this->setInteractive(false);
    }
    QGraphicsView::mousePressEvent(event);
}

void QGraphicsVolumeView::mouseReleaseEvent(QMouseEvent* event) {
    if (event->buttons() & Qt::RightButton) {
        this->setInteractive(true);
    }
    QGraphicsView::mouseReleaseEvent(event);
}
