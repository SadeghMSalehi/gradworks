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
#include <QtAlgorithms>
#include <QPolygonF>
#include "qtypedef.h"

using namespace pi;


QGraphicsVolumeView::QGraphicsVolumeView(QWidget* parent): QGraphicsView(parent) {
    setScene(&_scene);

    _thumbsWidth = 128;
    _columnCount = 10000;
    
    _xFlipped = false;
    _yFlipped = false;

//    _displayId = 0;
//    _displayReference = true;
//    _manualIntensityScaling = false;
    _directionCache = Unknown;

    _airImages = NULL;
    _currentSliceMarker = NULL;

    QGraphicsView::setDragMode(QGraphicsView::ScrollHandDrag);
//    setBackgroundBrush(QBrush(QPixmap::fromImage(QImage(QString::fromUtf8(":/Icons/Images/backgroundPattern.jpg")))));
}

void QGraphicsVolumeView::setImageFlip(bool xFlip, bool yFlip) {
    _xFlipped = xFlip;
    _yFlipped = yFlip;
}


void QGraphicsVolumeView::setDisplayCollection(pi::AIRDisplayCollection *images, bool useNavigationImage) {
    if (images == NULL) {
        return;
    }
    if (_airImages == images) {
        return;
    }
    clear();

    _airImages = images;
    _useNavigationImage = useNavigationImage;

    if (_airImages->GetReferenceId() >= 0) {
        setVolumeToShow(_airImages->GetReferenceId());
    }
}

void QGraphicsVolumeView::fitToImage(int sliceIdx, int volumeId) {
    if (_airImages != NULL && _airImages->Count() > 0) {
        if (volumeId < 0) {
            if (_volumeDisplays.contains(volumeId)) {
                QGraphicsPixmapItem* item = _volumeDisplays[volumeId].GetSliceData<QGraphicsPixmapItem>(sliceIdx);
                if (item != NULL) {
                    fitInView(item, Qt::KeepAspectRatio);
                }
            }
        } else {
            QRectF sceneRect = this->sceneRect();
            if (_volumeDisplays.size() > 0) {
                QGraphicsPixmapItem* item = _volumeDisplays[0].GetSliceData<QGraphicsPixmapItem>(sliceIdx);
                QRectF viewport = item->boundingRect();
                if (item != NULL) {
                    if (sceneRect.width() > sceneRect.height()) {
                        viewport.setHeight(sceneRect.height());
                    } else {
                        viewport.setWidth(sceneRect.width());
                    }
                }
            }
        }
    }
}

void QGraphicsVolumeView::clear() {
    _scene.clear();
//    _sliceCache.clear();
//    _slicePixmaps.clear();
//    _displayImages.clear();
    _workingSet.clear();
    _directionCache = Unknown;
//    _volumeCache = NULL;
//    _volumeSource = NULL;
    _currentSliceMarker = NULL;
    _volumeDisplays.clear();
}

void QGraphicsVolumeView::setVolumeToShow(int i) {
    if (i >= _airImages->Count() || i < 0) {
        return;
    }
    if (_volumeDisplays.contains(i)) {
        if (_volumeDisplays[i].Has(_airImages->at(i))) {
            return;
        } else {
            _volumeDisplays.remove(i);
        }
    }
    AIRVolumeDisplay newVolume;
    newVolume.SetDisplay(_airImages->at(i));

    _volumeDisplays[i] = newVolume;
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

//bool QGraphicsVolumeView::checkVolumeCache() {
//    AIRImageDisplay& refImg = _airImages->GetReference();
//    if (refImg.srcImg.IsNull()) {
//        _volumeCache = NULL;
//        _volumeSource = NULL;
//        return false;
//    }
//
//    // check if source is different pointer or the source itself has changed
//    if (_volumeSource != refImg.srcImg || _volumeSourceModifiedTime != refImg.srcImg->GetMTime()) {
//        clear();
//        _volumeSource = refImg.srcImg;
//        _volumeSourceModifiedTime = refImg.srcImg->GetMTime();
//    }
//
//    AIRImage::SpacingType spacing = refImg.srcSpacing;
//    AIRImage::SizeType size = refImg.srcImg->GetBufferedRegion().GetSize();
//
//    _rescaleFactor = 1;
//    if (_thumbsWidth != 0) {
//        _rescaleFactor = size[0] / _thumbsWidth;
//    }
//
//    if (std::abs(_rescaleFactor - 1) < 0.1) {
//        _volumeCache = refImg.srcImg;
//        return _volumeCache.IsNotNull();
//    }
//
//    fordim(k) {
//        spacing[k] = spacing[k] * _rescaleFactor;
//        size[k] = size[k] / _rescaleFactor;
//    }
//
//    typedef itk::ResampleImageFilter<AIRImage, AIRImage> ResampleFilter;
//    ResampleFilter::Pointer resampler = ResampleFilter::New();
//    resampler->SetInput(refImg.srcImg);
//    resampler->SetOutputParametersFromImage(refImg.srcImg.GetPointer());
//    resampler->SetOutputSpacing(spacing);
//    resampler->SetSize(size);
//    resampler->Update();
//    _volumeCache = resampler->GetOutput();
//    _volumeCache->DisconnectPipeline();
//
//    return _volumeCache.IsNotNull();
//}

//bool QGraphicsVolumeView::checkSliceCache() {
//    AIRImageDisplay& refImg = _airImages->GetReference();
//    SliceDirectionEnum dir = refImg.GetResampleDirection(0);
//
//    // fire only when the slice direction is changed
//    if (_directionCache == dir) {
//        return true;
//    }
//
//    _directionCache = dir;
//    return false;
//}

//bool QGraphicsVolumeView::updateSource() {
//    if (_directionCache == Unknown) {
//        return false;
//    }
//    AIRImage::RegionType region = _volumeCache->GetBufferedRegion();
//
//    _sliceCache.clear();
//    _sliceCache.reserve(region.GetSize(_directionCache));
//
//    typedef itk::ExtractImageFilter<AIRImage, AIRImage> ExtractFilter;
//
//    const int nSlices = region.GetSize(_directionCache);
//    for (int i = 0; i < nSlices; i++) {
//        ExtractFilter::Pointer extract = ExtractFilter::New();
//        extract->SetInput(_volumeCache);
//        region.SetIndex(_directionCache, i);
//        region.SetSize(_directionCache, 1);
//        extract->SetExtractionRegion(region);
//        extract->Update();
//        _sliceCache.push_back(extract->GetOutput());
//        _sliceCache.back()->DisconnectPipeline();
//    }
//    return true;
//}

void QGraphicsVolumeView::updateDisplay(int volumeId) {
    if (this->isHidden()) {
        return;
    }
    
    if (_airImages == NULL) {
        return;
    }

    _directionCache = _airImages->GetReference().GetResampleDirection(0);

    QIntList showingVolumes;
    if (volumeId < 0) {
        showingVolumes = _volumeDisplays.keys();
        qSort(showingVolumes);
    } else {
        if (!_volumeDisplays.contains(volumeId)) {
            return;
        }
        showingVolumes.append(volumeId);
    }

    _currentSliceMarker = NULL;

    QFont sliceIndexFont("Courier", 20);
    
    int volumeCount = 0;
    QIntList::ConstIterator iter;
    for (iter = showingVolumes.constBegin(); iter != showingVolumes.constEnd(); iter++) {
        int id = *iter;
        if (!_airImages->IsValidId(id)) {
            _volumeDisplays.remove(id);
            continue;
        }
        AIRImageDisplay* src = _airImages->at(id);
        if (!_volumeDisplays[id].Has(src)) {
            _volumeDisplays[id].SetDisplay(src);
        }
        if (_volumeDisplays[id].UpdateSlice(_directionCache, _useNavigationImage)) {
            const int w = _volumeDisplays[id].Width();
            const int h = _volumeDisplays[id].Height();
            const int s = _volumeDisplays[id].Count();

            for (int i = 0; i < s; i++) {
                int realSliceIdx = i;
                if (_useNavigationImage) {
                    realSliceIdx = _airImages->at(id)->GetNavigationImage().GetOriginalIndex(i);
                }

                int colPos = i * w;
                int rowPos = volumeCount * h;
                QGraphicsPixmapItem* item = _volumeDisplays[id].GetSliceData<QGraphicsPixmapItem>(i);
                QPixmap pixmap = QPixmap::fromImage(QImage(_volumeDisplays[id].GetColorImageBuffer(i), w, h, QImage::Format_ARGB32));
                if (item == NULL) {
                    item = _scene.addPixmap(QPixmap::fromImage(QImage(_volumeDisplays[id].GetColorImageBuffer(i), w, h, QImage::Format_ARGB32)));
                    _volumeDisplays[id].SetSliceData(i, item);
                    item->setPos(colPos, rowPos);
//                    item->setFlags(QGraphicsItem::ItemIsSelectable);
                    item->setData(SliceIndex, QVariant(i));
                    item->setData(AnnotationType, QVariant(SliceImage));
                    item->setData(RealSliceIndex, QVariant(realSliceIdx));
                } else {
                    item->setPixmap(pixmap);
                }

                QGraphicsTextItem* text = new QGraphicsTextItem(QString("#%1.%2").arg(id).arg(realSliceIdx), item);
                text->setFont(sliceIndexFont);
                text->setPos(3, 3);
                text->setZValue(1);
                text->setDefaultTextColor(Qt::yellow);
                
                if (_workingSet.contains(i)) {
                    addWorkingSetItem(item);
                }
            }
        }
        volumeCount++;
    }
}

void QGraphicsVolumeView::addWorkingSetItem(QGraphicsPixmapItem *sliceItem) {
    const static int markerSize = 9;
    QRectF rect = sliceItem->boundingRect();
    QRectF ellipse = QRectF(rect.width() - markerSize - 3, 3, markerSize, markerSize);
    QGraphicsEllipseItem* rectItem = new QGraphicsEllipseItem(ellipse, sliceItem);
    rectItem->setZValue(1);
    rectItem->setPen(QPen(Qt::green, 3));
    rectItem->setData(AnnotationType, QVariant(WorkingSet));
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
            removeWorkingSetItem(sliceIdx);
        } else {
            _workingSet.insert(sliceIdx);
            addWorkingSetItem(sliceItem);
        }
    }
    this->setInteractive(true);
}

void QGraphicsVolumeView::removeWorkingSetItem(int idx) {
    QList<QGraphicsItem*> children = _volumeDisplays[0].GetSliceData<QGraphicsPixmapItem>(idx)->childItems();
    QList<QGraphicsItem*>::ConstIterator iter = children.begin();
    for (; iter != children.end(); iter++) {
        if ((*iter)->data(AnnotationType).value<int>() == WorkingSet) {
            _scene.removeItem(*iter);
        }
    }
    _workingSet.remove(idx);
}

void QGraphicsVolumeView::clearWorkingSet() {
    QIntSet::ConstIterator iter = _workingSet.constBegin();
    for (; iter != _workingSet.constEnd(); iter++) {
        removeWorkingSetItem(*iter);
    }
    _workingSet.clear();
    this->setInteractive(true);
}

void QGraphicsVolumeView::currentSliceChanged(int slice) {
    if (this->isHidden()) {
        return;
    }
    if (_directionCache == Unknown) {
        return;
    }
    if (_rescaleFactor <= 0) {
        return;
    }

    int sliceIdx = slice;
    if (_useNavigationImage) {
        sliceIdx = _airImages->at(0)->GetNavigationImage().GetIndexFromOriginal(sliceIdx);
    }
    if (sliceIdx < 0 || sliceIdx >= _volumeDisplays[0].Count()) {
        return;
    }

    QGraphicsPixmapItem* currentSlicePixmap = _volumeDisplays[0].GetSliceData<QGraphicsPixmapItem>(sliceIdx);
    if (_currentSliceMarker == NULL) {
        _currentSliceMarker = new QGraphicsRectItem(currentSlicePixmap->boundingRect(), currentSlicePixmap);
        _currentSliceMarker->setPen(QPen(Qt::green, 3));
        _currentSliceMarker->setData(AnnotationType, SliceMarker);
    } else {
        _currentSliceMarker->setParentItem(currentSlicePixmap);
    }
}

void QGraphicsVolumeView::moveToVolume(int i) {
    if (_volumeDisplays.contains(i)) {
        QGraphicsPixmapItem* firstItem = _volumeDisplays[i].GetSliceData<QGraphicsPixmapItem>(0);
        if (firstItem == NULL) {
            return;
        }

        QRect contentsRect = viewport()->contentsRect();
        QPolygonF visibleScene = mapToScene(contentsRect);
        QRectF visibleRect = visibleScene.boundingRect();
        QPointF firstPosition = firstItem->pos();

        visibleRect.moveTop(firstPosition.y());
        this->ensureVisible(visibleRect, 0, 0);
//        QPoint displacement = mapFromScene(0, firstPosition.y() - visibleRect.top());
//        this->scroll(0, displacement.y());
    }
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

void QGraphicsVolumeView::mouseDoubleClickEvent(QMouseEvent* event) {
    if (event->buttons() & Qt::LeftButton) {
        QPointF pos = mapToScene(event->pos());
        QGraphicsItem* item = _scene.itemAt(pos);
        QGraphicsPixmapItem* sliceView = dynamic_cast<QGraphicsPixmapItem*>(item);
        if (sliceView != NULL) {
            int slice = sliceView->data(RealSliceIndex).value<int>();
            emit sliceDoubleClicked(slice);
        }
        return;
    }
    QGraphicsView::mouseDoubleClickEvent(event);

}


void QGraphicsVolumeView::dragEnterEvent(QDragEnterEvent *event) {
    if (event->mimeData()->hasUrls()) {
        UrlList urls = event->mimeData()->urls();
        UrlList::ConstIterator iter;
        for (iter = urls.constBegin(); iter != urls.constEnd(); iter++) {
            const QUrl& url = (*iter);
            if (url.scheme() != "file") {
                event->setAccepted(false);
                return;
            }
        }
        event->setAccepted(true);
        event->acceptProposedAction();
    }
}

void QGraphicsVolumeView::dropEvent(QDropEvent *event) {
    QList<QString> fileNames;
    if (event->mimeData()->hasUrls()) {
        UrlList urls = event->mimeData()->urls();
        const QUrl& url = urls[0];
        if (url.scheme() != "file") {
            return;
        }
        event->acceptProposedAction();
        QString filePath = url.path();
        emit fileDropped(filePath);
    }
}