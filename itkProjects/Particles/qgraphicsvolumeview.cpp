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
    _scene = new QGraphicsScene();
    setScene(_scene);

    _thumbsWidth = 128;
    _columnCount = 10000;
    
    _xFlipped = false;
    _yFlipped = false;
    _showAll = true;
    _directionCache = IJ;

    _airImages = NULL;
    _currentSliceMarker = NULL;

    QGraphicsView::setDragMode(QGraphicsView::ScrollHandDrag);
}


void QGraphicsVolumeView::flipLR(bool toggle) {
    _xFlipped = toggle;
    updatePixmaps();
}


void QGraphicsVolumeView::flipUD(bool toggle) {
    _yFlipped = toggle;
    updatePixmaps();
}


void QGraphicsVolumeView::directionChanged(pi::SliceDirectionEnum dir) {
    if (_directionCache == dir) {
        return;
    }

    _directionCache = dir;
    delete _scene;

    _scene = new QGraphicsScene();
    setScene(_scene);

    updateDisplay();
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

    if (!_airImages->IsEmpty()) {
    }
}

void QGraphicsVolumeView::fitToImage(int sliceIdx, int volumeId) {
    if (_airImages != NULL && _airImages->Count() > 0) {
        // clear scene rect to recompute bounding box
        QRectF sceneRect = scene()->sceneRect();
        if (_volumeDisplays.size() > 0) {
            QGraphicsItem* item = _volumeDisplays[0].GetSliceData<QGraphicsItem>(sliceIdx);
            QRectF viewport = item->boundingRect();
            if (item != NULL) {
                if (sceneRect.width() > sceneRect.height()) {
                    viewport.setHeight(sceneRect.height());
                } else {
                    viewport.setWidth(sceneRect.width());
                }
            }
            fitInView(viewport, Qt::KeepAspectRatio);
            centerOn(item->pos().x() + viewport.width()/2, item->pos().y() + viewport.height()/2);
        }
    }
}

void QGraphicsVolumeView::clear() {
    _workingSet.clear();
    _directionCache = IJ;
    _currentSliceMarker = NULL;
    _volumeDisplays.clear();
    delete _scene;
    
    _scene = new QGraphicsScene();
    setScene(_scene);
}

/*
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
*/

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


void QGraphicsVolumeView::updateDisplay(int volumeId) {
    if (this->isHidden()) {
        return;
    }
    
    if (_airImages == NULL) {
        return;
    }

    _currentSliceMarker = NULL;

    QFont sliceIndexFont("Courier", 20);

    int showingVolumeCount = _showAll ? _airImages->Count() : 1;
    for (int id = 0; id < showingVolumeCount; id++) {
        AIRImageDisplay src = _airImages->at(id);
        if (_volumeDisplays.size() <= id) {
            _volumeDisplays.push_back(AIRVolumeDisplay());
        }
        if (!_volumeDisplays[id].Has(src)) {
            _volumeDisplays[id].SetDisplay(src);
        }
        if (_volumeDisplays[id].UpdateSlice(_directionCache, _useNavigationImage)) {
            const int w = _volumeDisplays[id].Width();
            const int h = _volumeDisplays[id].Height();
            const int s = _volumeDisplays[id].Count();

            for (int i = 0; i < s; i++) {
                int realSliceIdx = i;
                int colPos = i * w;
                int rowPos = id * h;

                QGraphicsRectItem* item = _volumeDisplays[id].GetSliceData<QGraphicsRectItem>(i);
                QGraphicsPixmapItem* pixmapItem = NULL;
                uchar* colorPointer = _volumeDisplays[id].GetColorImageBuffer(i);

                QPixmap pixmap = QPixmap::fromImage(QImage(colorPointer, w, h, QImage::Format_ARGB32));
                if (item == NULL) {
                    item = new QGraphicsRectItem(QRect(0, 0, w, h));
                    item->setPen(Qt::NoPen);
                    item->setBrush(Qt::NoBrush);
                    item->setPos(colPos, rowPos);
                    item->setData(SliceIndex, QVariant(i));
                    item->setData(AnnotationType, QVariant(SliceImage));
                    item->setData(RealSliceIndex, QVariant(realSliceIdx));

                    pixmapItem = new QGraphicsPixmapItem(QPixmap::fromImage(QImage(_volumeDisplays[id].GetColorImageBuffer(i), w, h, QImage::Format_ARGB32)), item);
                    pixmapItem->setZValue(1);

                    QGraphicsTextItem* text = new QGraphicsTextItem(QString("#%1.%2").arg(id).arg(realSliceIdx), item);
                    text->setFont(sliceIndexFont);
                    text->setPos(3, 3);
                    text->setZValue(2);
                    text->setDefaultTextColor(Qt::yellow);

                    _volumeDisplays[id].SetSliceData(i, item);
                    scene()->addItem(item);
                } else {
                    QGraphicsPixmapItem* pixmapItem = (QGraphicsPixmapItem*) item->childItems()[0];
                    pixmapItem->setPixmap(pixmap);
                }

                if (_workingSet.contains(i)) {
                    addWorkingSetItem(item);
                }
            }
        }
    }
    updatePixmaps();
}

void QGraphicsVolumeView::updatePixmaps() {
    for (int i = 0; i < _volumeDisplays.size(); i++) {
        const int w = _volumeDisplays[i].Width();
        const int h = _volumeDisplays[i].Height();
        const int s = _volumeDisplays[i].Count();

        for (int j = 0; j < s; j++) {
            AIRVolumeDisplay& display = _volumeDisplays[i];
            QGraphicsRectItem* item = display.GetSliceData<QGraphicsRectItem>(j);
            if (item != NULL) {
                QList<QGraphicsItem*> childItems = item->childItems();
                if (childItems.size() > 0) {
                    QGraphicsPixmapItem* pixmapItem = (QGraphicsPixmapItem*) childItems[0];
                    QTransform pixmapTransform;
                    if (_xFlipped && !_yFlipped) {
                        pixmapTransform.setMatrix(-1, 0, 0, 0, 1, 0, w, 0, 1);
                        pixmapItem->setTransform(pixmapTransform);
                    } else if (!_xFlipped && _yFlipped) {
                        pixmapTransform.setMatrix(1, 0, 0, 0, -1, 0, 0, h, 1);
                        pixmapItem->setTransform(pixmapTransform);
                    } else if (_xFlipped && _yFlipped) {
                        pixmapTransform.setMatrix(-1, 0, 0, 0, -1, 0, w, h, 1);
                        pixmapItem->setTransform(pixmapTransform);
                    } else {
                        pixmapItem->resetTransform();
                    }

                }
            }
        }
    }
}

void QGraphicsVolumeView::addWorkingSetItem(QGraphicsItem *sliceItem) {
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
    QList<QGraphicsItem*> selections = scene()->selectedItems();
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
            scene()->removeItem(*iter);
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
        //sliceIdx = _airImages->at(0)->GetNavigationImage().GetIndexFromOriginal(sliceIdx);
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
    if (i >= _volumeDisplays.size()) {
        return;
    }

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
        QGraphicsItem* item = scene()->itemAt(pos);
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