//
//  QGraphicsGuideView.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/24/13.
//
//

#include "qgraphicsguideview.h"
#include "QMouseEvent"
#include "QGraphicsScene"
#include <QPainterPath>
#include <QGraphicsPathItem>
#include <QPainter>
#include "qgraphicscompositeimageitem.h"
#include "piImageIO.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "piDelegateImageFilter.h"

using namespace std;
using namespace pi;


uint __labelColors[254] = { 0,
    0xff930f2b, 0xfff4da0b, 0xff1425b3, 0xffb5e489, 0xffc33bf4, 0xffdb1d6b, 0xff847a79, 0xffe0763c, 0xff78eaa7,
    0xffb0c74b, 0xff3b74a8, 0xfff6ebe8, 0xff3b7906, 0xff2e0cd0, 0xff3ffadc, 0xfff06b12, 0xff174909, 0xffbf9121,
    0xff9646af, 0xff2be490, 0xffae0611, 0xff754d51, 0xffe27c01, 0xff3a9708, 0xffc6853c, 0xfff84707, 0xff73dda6,
    0xff8392b4, 0xff54f584
};
uint __tmpLabelColors[254] = { 0 };
uint __maxLabelColors = 30;
uint* __labelColorsPtr = __labelColors;

#pragma mark -
#pragma mark Auxilirary Class

template <class T>
class LabelPropagation {
public:
    typedef typename T::Pointer TPointer;
    
    TPointer& _input;
    TPointer& _output;
    IntVector& _target;
    
    LabelPropagation(AIRLabel::Pointer& input, AIRLabel::Pointer& output, IntVector& target): _input(input), _output(output), _target(target) {
        
    }

    void Process() {
        IntVector::const_iterator iter(_target.begin());
        for (; iter != _target.end(); iter++) {
            Process(*iter);
        }
    }
    
    void Process(int sliceIdx) {
        AIRLabel::RegionType region = _input->GetBufferedRegion();
        for (int i = 0; i < DIMENSIONS; i++) {
            if (region.GetSize(i) == 1) {
                region.SetIndex(i, sliceIdx);
            }
        }
        typedef typename itk::ImageRegionIterator<T> IterType;
        IterType srcIter(_input, _input->GetBufferedRegion());
        IterType dstIter(_output, region);
        for (srcIter.GoToBegin(), dstIter.GoToBegin(); !srcIter.IsAtEnd() && !dstIter.IsAtEnd(); ++srcIter, ++dstIter) {
            dstIter.Set(srcIter.Get());
        }
    }
};


#pragma mark -
#pragma mark Main Implementation

QGraphicsGuideView::QGraphicsGuideView(QWidget* parent): QGraphicsView(parent) {
    _currentMode = TransformMode;
    _labelItem = NULL;
    _brushCursor = NULL;
    _brushRadius = 10.5;

    _volumeSliceIdx = -1;
    _volumeSliceDirection = Unknown;

    _foregroundLabel = 1;
    _backgroundLabel = 255;
}

QGraphicsGuideView::~QGraphicsGuideView() {
}

void QGraphicsGuideView::drawingMode(ToolType toolType) {
    _currentMode = DrawingMode;
    _currentTool = toolType;
    updateLabelItem();
}

void QGraphicsGuideView::transformMode() {
    _currentMode = TransformMode;
}

AIRLabelSlice QGraphicsGuideView::getLabelSlice() {
    return _volumeSlice;
}

void QGraphicsGuideView::mousePressEvent(QMouseEvent* event) {
    if (event->buttons() & Qt::LeftButton) {
        if (DrawingMode == _currentMode) {
            if (_labelImage.size().width() == 0) {
                segmentationCleared();
            }
            event->accept();
            switch (_currentTool) {
                case FreePath:
                    pathStartEvent(event);
                    setCursor(Qt::CrossCursor);
                    return;
                case FreeBrush:
                case EraseBrush:
                    setCursor(Qt::BlankCursor);
                    brushStartEvent(event);
                    return;
                default:
                    break;
            }
        }
    }
    QGraphicsView::mousePressEvent(event);
}

void QGraphicsGuideView::mouseMoveEvent(QMouseEvent* event) {
    if (DrawingMode == _currentMode) {
        event->accept();
        switch (_currentTool) {
            case FreePath:
                pathMoveEvent(event);
                return;
            case FreeBrush:
            case EraseBrush:
                brushMoveEvent(event);
                return;
            default:
                break;
        }
    }

    QGraphicsView::mouseMoveEvent(event);
}

void QGraphicsGuideView::mouseReleaseEvent(QMouseEvent* event) {
    if (DrawingMode == _currentMode) {
        event->accept();
        switch (_currentTool) {
            case FreePath:
                pathEndEvent(event);
                unsetCursor();
                return;
            case FreeBrush:
            case EraseBrush:
                brushEndEvent(event);
                unsetCursor();                
                return;
            default:
                break;
        }
    }

    QGraphicsView::mouseReleaseEvent(event);
}

void QGraphicsGuideView::pathStartEvent(QMouseEvent *event) {
    QPointF pos = mapToScene(event->pos());
    _drawingPressPoint = pos;
    _drawingPath = QPainterPath(_drawingPressPoint);

    if (this->scene() != NULL) {
        _drawingPathItem = this->scene()->addPath(_drawingPath);
        _drawingPathItem->setOpacity(0.2);
    }
    _drawingPathItem->setPen(QPen(Qt::green,1));
}

void QGraphicsGuideView::pathMoveEvent(QMouseEvent *event) {
    QPointF pos = mapToScene(event->pos());
    _drawingPath.lineTo(pos);
    _drawingPathItem->setPath(_drawingPath);
    _drawingPathItem->update();
}

void QGraphicsGuideView::pathEndEvent(QMouseEvent *event) {
    QPainter painter(&_drawingCanvas);
    painter.setBrush(QBrush(Qt::yellow,Qt::SolidPattern));
    painter.setPen(Qt::NoPen);
    painter.drawPath(_drawingPath);
    painter.end();

    this->scene()->removeItem(_drawingPathItem);

    canvasToSlice();
    sliceToLabelImage();
    updateLabelItem();
    sliceToVolume();
}

void QGraphicsGuideView::brushStartEvent(QMouseEvent *event) {

    QPoint topLeft(event->pos().x()-_brushRadius,event->pos().y()-_brushRadius);
    QPoint bottomRight(event->pos().x()+_brushRadius,event->pos().y()+_brushRadius);
    QPointF sceneTopLeft = mapToScene(topLeft);
    QPointF sceneBottomRight = mapToScene(bottomRight);
    QRectF sceneBrushRect(sceneTopLeft, sceneBottomRight);

    if (_currentTool == FreeBrush) {
        _brushCursor = this->scene()->addEllipse(0,0,sceneBrushRect.width(), sceneBrushRect.height());
    } else if (_currentTool == EraseBrush ) {
        _brushCursor = this->scene()->addEllipse(0,0,sceneBrushRect.width(), sceneBrushRect.height());
    }

    _brushCursor->setPos(sceneTopLeft);
    _brushCursor->setPen(QPen(Qt::white, 1));
    _brushCursor->setZValue(10);
    _brushCursor->setOpacity(0.5);

    paintBrushEllipse(sceneBrushRect);
    canvasToSlice();
    sliceToLabelImage();
    updateLabelItem();
}

void QGraphicsGuideView::brushMoveEvent(QMouseEvent *event) {
    if (_brushCursor == NULL) {
        return;
    }
    QPoint topLeft(event->pos().x()-_brushRadius,event->pos().y()-_brushRadius);
    QPoint bottomRight(event->pos().x()+_brushRadius,event->pos().y()+_brushRadius);
    QPointF sceneTopLeft = mapToScene(topLeft);
    QPointF sceneBottomRight = mapToScene(bottomRight);
    QRectF sceneBrushRect(sceneTopLeft, sceneBottomRight);

    _brushCursor->setPos(sceneTopLeft);

    paintBrushEllipse(sceneBrushRect);
    canvasToSlice();
    sliceToLabelImage();
    updateLabelItem();
}

void QGraphicsGuideView::brushEndEvent(QMouseEvent *event) {
    if (_brushCursor == NULL) {
        return;
    }
    sliceToVolume();
    this->scene()->removeItem(_brushCursor);
    _brushCursor = NULL;
}


#pragma mark -
#pragma mark Private Methods

void QGraphicsGuideView::paintBrushEllipse(QRectF& sceneBrushRect) {
    QPainter painter(&_drawingCanvas);
    painter.setBrush(QBrush(Qt::yellow,Qt::SolidPattern));
    painter.setPen(Qt::NoPen);
    painter.drawEllipse(sceneBrushRect);
    painter.end();
}

void QGraphicsGuideView::canvasToSlice() {
    if (_volumeSlice.IsNull()) {
        return;
    }

    const int h = _drawingCanvas.height();
    const int w = _drawingCanvas.width();

    uint color = 0;
    if (_currentTool != EraseBrush) {
        color = qRgba(255, 255, 0, 255);
    }

    QImage sliceImage(_volumeSlice.GetBufferPointer(), w, h, QImage::Format_Indexed8);

    if (_currentTool == EraseBrush) {
        ushort* src = (ushort*) _drawingCanvas.scanLine(0);
        uchar* dst = (uchar*) sliceImage.scanLine(0);

        for (int i = 0; i < w*h; i++) {
            if (*src > 0) {
                *src = *dst = 0;
            }
            src++;
            dst++;
        }
    } else {
        ushort* src = (ushort*) _drawingCanvas.scanLine(0);
        uchar* dst = (uchar*) sliceImage.scanLine(0);
        for (int i = 0; i < w*h; i++) {
            if (*src > 0 && (_backgroundLabel == 255 || *dst == _backgroundLabel)) {
                *dst = _foregroundLabel;
                *src = 0;
            }
            src++;
            dst++;
        }
    }
}



void QGraphicsGuideView::updateLabelItem() {
    if (this->scene() == NULL) {
        return;
    }
    if (_labelItem == NULL) {
        _labelItem = this->scene()->addPixmap(QPixmap::fromImage(_labelImage));
        _labelItem->setOpacity(0.5);
        _labelItem->setZValue(10);
    } else {
        _labelItem->setPixmap(QPixmap::fromImage(_labelImage));
        _labelItem->update();
    }
}

bool QGraphicsGuideView::isLabelLoaded() {
    return _labelVolume.IsNotNull();
}

void QGraphicsGuideView::sliceToVolume() {
    typedef itk::ImageRegionConstIteratorWithIndex<AIRLabel> IterType;
    IterType iter(_volumeSlice, _volumeSlice.GetRegion());
    for (iter.GoToBegin(); !iter.IsAtEnd(); ++iter) {
        _labelVolume->SetPixel(iter.GetIndex(), iter.Get());
    }
}

void QGraphicsGuideView::volumeToSlice() {
    if (_volumeSliceIdx < 0) {
        return;
    }
    _volumeSlice = ExtractSlice<AIRLabel>(_labelVolume, _volumeSliceIdx, _volumeSliceDirection);
}

void QGraphicsGuideView::sliceToLabelImage() {
    if (_volumeSlice.IsNull()) {
        return;
    }

    AIRClass* slicePixels = _volumeSlice.GetBufferPointer();

    const int w = _labelImage.width();
    const int h = _labelImage.height();

    for (int i = 0; i < h; i++) {
        uint* labelPixels = (uint*) _labelImage.scanLine(i);
        for (int j = 0; j < w; j++) {
            (*labelPixels) = __labelColorsPtr[*slicePixels];
            labelPixels++;
            slicePixels++;
        }
    }
}

#pragma mark -
#pragma mark Public Slots

void QGraphicsGuideView::segmentationCleared() {
    if (_volumeSlice.IsNull()) {
        return;
    }
    _volumeSlice.FillBuffer(0);
    sliceToLabelImage();
    sliceToVolume();
    updateLabelItem();
}

void QGraphicsGuideView::labelOpacityChanged(int n) {
    if (_labelItem != NULL) {
        _labelItem->setOpacity(n/255.0);
    }

}

void QGraphicsGuideView::sliceChanged(SliceDirectionEnum dir, int n) {
    if (dir == Unknown) {
        return;
    }

    _volumeSliceIdx = n;
    bool updateCanvas = false;
    if (_volumeSliceDirection != dir) {
        _volumeSliceDirection = dir;
        updateCanvas = true;
    }

    volumeToSlice();
    if (_volumeSlice.IsNull()) {
        return;
    }
    if (updateCanvas) {
        // when slice direction is changed
        SliceDisplay<AIRLabel> grid(_volumeSlice);
        _canvasSize = QSize(grid.Width(), grid.Height());
        _labelImage = QImage(_canvasSize, QImage::Format_ARGB32_Premultiplied);
        _labelImage.fill(0);
        _drawingCanvas = QImage(_canvasSize, QImage::Format_RGB16);
        _drawingCanvas.fill(0);

    }
    sliceToLabelImage();
    updateLabelItem();
}

void QGraphicsGuideView::labelVolumeChanged(AIRLabel::Pointer labelVolume) {
    _labelVolume = labelVolume;
    volumeToSlice();
}

void QGraphicsGuideView::createLabelVolumeIfNecessary(AIRImage::Pointer srcImg) {
    if (_labelVolume.IsNotNull()) {
        if (_labelVolume->GetBufferedRegion() == srcImg->GetBufferedRegion()) {
            // no need to change buffer
            // do not clean
            cout << "keep label volume in spite of new source volume:" << __FILE__ << endl;
            return;
        }
    }
    AIRLabel::Pointer labelVolume = __airImageIO.CastImageToS<AIRLabel>(srcImg);
    labelVolume->FillBuffer(0);
    labelVolumeChanged(labelVolume);
}

void QGraphicsGuideView::loadLabelVolume(QString& fileName, AIRImage::Pointer srcImg) {
    AIRLabel::Pointer labelVolume = __airLabelIO.ReadImage(fileName.toUtf8().data());
    if (_labelVolume.IsNull()) {
        labelVolumeChanged(labelVolume);
    } else {
        if (srcImg->GetBufferedRegion() != labelVolume->GetBufferedRegion()) {
            return;
        }
        labelVolumeChanged(labelVolume);
        volumeToSlice();
        sliceToLabelImage();
        updateLabelItem();
    }
}

void QGraphicsGuideView::saveLabelVolume(QString& fileName) {
    __airLabelIO.WriteImage(fileName.toUtf8().data(), _labelVolume);
}

void QGraphicsGuideView::propagateLabel(IntVector &targetSlices) {
    AIRLabel::Pointer label = _volumeSlice.GetImage();
    LabelPropagation<AIRLabel> algo(label, _labelVolume, targetSlices);
    algo.Process();
 }

void QGraphicsGuideView::setInteractionLabels(int f, int b) {
    _foregroundLabel = f;
    _backgroundLabel = b;
}

void QGraphicsGuideView::copyLabel() {
    if (_volumeSlice.IsNotNull()) {
        _copiedSlice = __airLabelIO.CopyImage(_volumeSlice.GetImage());
        cout << "copyLabel()" << endl;
    }
}

void QGraphicsGuideView::pasteLabel() {
    if (_copiedSlice.IsNotNull()) {
        AIRLabelSlice copied(_copiedSlice);
        if (_volumeSlice.IsNotNull()
            && _volumeSlice.Direction() == copied.Direction()
            && _volumeSlice.Width() == copied.Width()
            && _volumeSlice.Height() == copied.Height()) {

            AIRClass* src = _copiedSlice->GetBufferPointer();
            AIRClass* dst = _volumeSlice.GetImage()->GetBufferPointer();
            memcpy(dst, src, _copiedSlice->GetPixelContainer()->Size());
            sliceToVolume();
            sliceToLabelImage();
            updateLabelItem();
        } else {
            _copiedSlice = NULL;
        }
    }
}