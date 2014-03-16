//
//  QGraphicsImageItem.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/7/13.
//
//

#ifndef __ParticleGuidedRegistration__QGraphicsImageItem__
#define __ParticleGuidedRegistration__QGraphicsImageItem__

#include <iostream>
#include <itkImage.h>
#include <QGraphicsItem>
#include <QGraphicsSceneMouseEvent>
#include <QGraphicsSceneHoverEvent>
#include <QPainter>
#include <QWidget>

#include "piImageHistogram.h"


extern unsigned int __grays[256];
extern unsigned int __red2blue[256];
extern unsigned int __blue2red[256];


template <class T>
class QGraphicsItemInteraction {
public:
    virtual void mousePressed(T*,QGraphicsSceneMouseEvent*) = 0;
    virtual void mouseMoved(T*,QGraphicsSceneMouseEvent*) = 0;
    virtual void mouseReleased(T*,QGraphicsSceneMouseEvent*) = 0;
    virtual void hoverEntered(T*,QGraphicsSceneHoverEvent*) = 0;
    virtual void hoverMoved(T*,QGraphicsSceneHoverEvent*) = 0;
    virtual void hoverLeft(T*,QGraphicsSceneHoverEvent*) = 0;
};

template <class T>
class QGraphicsImageItem: public QGraphicsPixmapItem {
public:
    typedef T PixelType;
    typedef QVector<QPointF> GridCoord;
    typedef QGraphicsItemInteraction< QGraphicsImageItem<T> > InteractionType;

    enum Flags { NoFlip, LeftRight, UpDown };
    enum Mode { NoMode, GridMode };

    QGraphicsImageItem(QGraphicsItem* parent = NULL): QGraphicsPixmapItem(parent) {
        _range[0] = 0;
        _range[1] = 0;
        _showGrid = false;
        _interaction = NULL;
        setAcceptedMouseButtons(Qt::LeftButton);
    }
    virtual ~QGraphicsImageItem() {}
    virtual void paint(QPainter * painter, const QStyleOptionGraphicsItem * option, QWidget * widget = 0);

    T getRange(int i);
    void setRange(T min, T max);
    void computeRange();
    void setFlip(Flags flags);
    void showGrid(bool);
    void setInteraction(InteractionType* interaction);

    inline GridCoord& userGrids() { return _userGrids; }

    void refresh();
    
    void setImage(T* inputBuffer, int w, int h);

    // templated member function
    template <class S> void setImage(typename S::Pointer image, bool computeRange = false);
    template <class S> void generateUserGrids(typename S::Pointer transform);

    // static function
    static void convertToIndexed(T*, double, double, QImage&);


protected:
    void mousePressEvent(QGraphicsSceneMouseEvent* event);
    void mouseReleaseEvent(QGraphicsSceneMouseEvent* event);
    void mouseMoveEvent(QGraphicsSceneMouseEvent* event);
    void hoverEnterEvent(QGraphicsSceneHoverEvent* event);
    void hoverMoveEvent(QGraphicsSceneHoverEvent* event);
    void hoverLeaveEvent(QGraphicsSceneHoverEvent* event);

    void drawRegularGrid(QPainter*);
    void drawUserGrid(QPainter*);

private:
    T* _inputBuffer;
    T _range[2];
    bool _showGrid;
    QImage _grayImage;
    InteractionType* _interaction;
    GridCoord _userGrids;
};


template <class T>
void QGraphicsImageItem<T>::showGrid(bool onoff) {
    _showGrid = onoff;
    update();
}

template <class T>
void QGraphicsImageItem<T>::setInteraction(InteractionType *interaction) {
    _interaction = interaction;
    setAcceptHoverEvents(true);
}


template <class T>
void QGraphicsImageItem<T>::drawRegularGrid(QPainter* painter) {
    const int w = _grayImage.width();
    const int h = _grayImage.height();
    QPen pen = painter->pen();
    static QPen blackPen(Qt::yellow, 1), dashPen(Qt::yellow, 1, Qt::DotLine);
    blackPen.setCosmetic(true);
    dashPen.setCosmetic(true);
    painter->setPen(blackPen);

    const int hSpacing = h / 20.0;
    int n = 0;
    for (int i = 0; i < h; i += hSpacing, n++) {
        if (n % 5 == 0) {
            painter->setPen(blackPen);
        } else {
            painter->setPen(dashPen);
        }
        painter->drawLine(0, i, w-1, i);
    }
    painter->setPen(blackPen);
    painter->drawLine(0, h-1, w, h-1);

    n = 0;
    for (int i = 0; i < w; i += hSpacing, n++) {
        if (n % 5 == 0) {
            painter->setPen(blackPen);
        } else {
            painter->setPen(dashPen);
        }
        painter->drawLine(i, 0, i, h-1);
    }
    painter->setPen(blackPen);
    painter->drawLine(w-1, 0, w-1, h);
}

template <class T>
void QGraphicsImageItem<T>::drawUserGrid(QPainter *painter) {
    const int w = _grayImage.width();
    const int h = _grayImage.height();

    if (_userGrids.size() != w * h) {
        std::cout << __FILE__ << ": the size of user grid doesn't match with the image" << std::endl;
        return;
    }

    static QPen blackPen(Qt::yellow, 1), dashPen(Qt::yellow, 1, Qt::DotLine);
    blackPen.setCosmetic(true);
    dashPen.setCosmetic(true);

    const int spacing = h / 20;
    int wcnt = 0;
    for (int i = 0; i < w; i+= spacing, wcnt++) {
        int hcnt = 0;
        for (int j = 0; j < h; j+=spacing, hcnt++) {
            QPointF& top = _userGrids[j * w + i];

            int rightIdx = i + spacing;
            if (rightIdx >= w) {
                rightIdx = w - 1;
            }
            QPointF& right = _userGrids[j * w + rightIdx];

            int bottomIdx = (j + spacing);
            if (bottomIdx >= h) {
                bottomIdx = h - 1;
            }
            QPointF& bottom = _userGrids[w * bottomIdx + i];

            if (hcnt % 5 == 0) {
                painter->setPen(blackPen);
            } else {
                painter->setPen(dashPen);
            }
            painter->drawLine(top, right);

            if (wcnt % 5 == 0) {
                painter->setPen(blackPen);
            } else {
                painter->setPen(dashPen);
            }
            painter->drawLine(top, bottom);
        }
    }

    painter->setPen(blackPen);
    for (int i = 0; i < w - 1; i++) {
        QPointF& top = _userGrids[(h - 1) * w + i];
        QPointF& right = _userGrids[(h - 1) * w + i + 1];
        painter->drawLine(top, right);
    }

    for (int j = 1; j < h - 1; j++) {
        QPointF& top = _userGrids[(j - 1) * w + w - 1];
        QPointF& bottom = _userGrids[j * w + w - 1];
        painter->drawLine(top, bottom);
    }
}

template <class T>
void QGraphicsImageItem<T>::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget* widget) {
    QGraphicsPixmapItem::paint(painter, option, widget);
    if (_showGrid) {
        if (_userGrids.size() == 0) {
            drawRegularGrid(painter);
        } else {
            drawUserGrid(painter);
        }
    }
}

template <class T>
void QGraphicsImageItem<T>::setImage(T* inputBuffer, int w, int h) {
    _inputBuffer = inputBuffer;
    if (w != _grayImage.width() || h != _grayImage.height()) {
        _grayImage = QImage(w, h, QImage::Format_Indexed8);
        _grayImage.setColorCount(256);
        for (int i = 0; i < 256; i++) {
            _grayImage.setColor(i, __grays[i]);
        }
    }
}


template <class T> template <class S>
void QGraphicsImageItem<T>::setImage(typename S::Pointer image, bool computeRange) {
    if (image.IsNull()) {
        std::cout << __FILE__ << ": image is null!" << std::endl;
        return;
    }
    typename S::SizeType sz = image->GetBufferedRegion().GetSize();
    int w = 0;
    int h = 0;
    if (S::ImageDimension == 2) {
        w = sz[0];
        h = sz[1];
    } else {
        for (int i = 0; i < S::ImageDimension; i++) {
            if (sz[i] > 1 && w == 0) {
                w = sz[i];
            } else {
                h = sz[i];
                break;
            }
        }
    }
    if (w == 0 && h == 0) {
        return;
    }
    if (computeRange) {
        const int nelems = w * h;
        T* buff = image->GetBufferPointer();
        int i = 0;
        for (i = 0; i < nelems; i++, buff++) {
            if (*buff != 0) {
                _range[0] = _range[1] = *buff;
                buff++;
                break;
            }
        }
        // FIXME : test if it works well
        i = 0;
        for (; i < nelems; i++, buff++) {
            if (*buff == 0) {
                continue;
            }
            _range[0] = std::min(*buff, _range[0]);
            _range[1] = std::max(*buff, _range[1]);
        }

        std::cout << "Range: " << _range[0] << "," << _range[1] << std::endl;
    }
    setImage(image->GetBufferPointer(), w, h);
}

template <class T> template <class S>
void QGraphicsImageItem<T>::generateUserGrids(typename S::Pointer transform) {
    _userGrids.clear();

    const int w = _grayImage.width();
    const int h = _grayImage.height();

    for (int j = 0; j < h; j++) {
        for (int i = 0; i < w; i++) {
            typename S::InputPointType pIn;
            pIn[0] = i;
            pIn[1] = j;
            typename S::OutputPointType pOut = transform->TransformPoint(pIn);
            _userGrids.append(QPointF(pOut[0], pOut[1]));
        }
    }
}

template <class T>
T QGraphicsImageItem<T>::getRange(int i) {
    return _range[i];
}


template <class T>
void QGraphicsImageItem<T>::setRange(T min, T max) {
    _range[0] = min;
    _range[1] = max;
}

template <class T>
void QGraphicsImageItem<T>::computeRange() {
    const int w = _grayImage.width();
    const int h = _grayImage.height();
    T* buff = _inputBuffer;
    _range[0] = *buff;
    _range[1] = *buff;
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            _range[0] = std::min(_range[0], *buff);
            _range[1] = std::max(_range[1], *buff);
            buff++;
        }
    }
}

template <class T>
void QGraphicsImageItem<T>::setFlip(Flags flags) {
    if (flags == LeftRight) {
        QTransform transform;
        transform.setMatrix(-1, 0, 0, 0, 1, 0, _grayImage.width(), 0, 1);
        this->setTransform(transform);
    } else if (flags == UpDown) {
        QTransform transform;
        transform.setMatrix(1, 0, 0, 0, -1, 0, 0, _grayImage.height(), 1);
        this->setTransform(transform);
    } else {
        QTransform transform;
        transform.reset();
        this->setTransform(transform);
    }
}

template <class T>
void QGraphicsImageItem<T>::refresh() {
    if (_inputBuffer == NULL) {
        return;
    }

    const T pixelRange = _range[1] - _range[0];
    const T pixelMin = _range[0];

    convertToIndexed(_inputBuffer, pixelMin, pixelRange, _grayImage);

    QPixmap pixmap = QPixmap::fromImage(_grayImage);
    setPixmap(pixmap);
}

template <class T>
void QGraphicsImageItem<T>::convertToIndexed(T* inputBuffer,
                                               const double min, const double range, QImage& outputImage) {
    const int width = outputImage.width();
    for (int i = 0; i < outputImage.height(); i++) {
        uchar* outputBuffer = outputImage.scanLine(i);
        for (int j = 0; j < width; j++, outputBuffer++, inputBuffer++) {
            if (*inputBuffer > min + range) {
                *outputBuffer = 255u;
            } else {
                *outputBuffer = uchar(255u * (*inputBuffer-min)/(range));
            }
        }
    }
}

template <class T>
void QGraphicsImageItem<T>::mousePressEvent(QGraphicsSceneMouseEvent *event) {
    QPoint pos = event->pos().toPoint();
    int idx = pos.y() * _grayImage.width() + pos.x();
    std::cout << pos.x() << ", " << pos.y() << "[" << idx << "]: " << _inputBuffer[idx] << std::endl;

    if (_interaction) {
        _interaction->mousePressed(this, event);
        return;
    }

    return QGraphicsPixmapItem::mousePressEvent(event);
}

template <class T>
void QGraphicsImageItem<T>::mouseMoveEvent(QGraphicsSceneMouseEvent *event) {
    if (_interaction) {
        _interaction->mouseMoved(this, event);
        return;
    }
    return QGraphicsPixmapItem::mouseMoveEvent(event);
}

template <class T>
void QGraphicsImageItem<T>::mouseReleaseEvent(QGraphicsSceneMouseEvent *event) {
    if (_interaction) {
        _interaction->mouseReleased(this, event);
        return;
    }

    return QGraphicsPixmapItem::mouseReleaseEvent(event);
}


template <class T>
void QGraphicsImageItem<T>::hoverEnterEvent(QGraphicsSceneHoverEvent *event) {
    if (_interaction) {
        _interaction->hoverEntered(this, event);
        return;
    }
    return QGraphicsPixmapItem::hoverEnterEvent(event);
}

template <class T>
void QGraphicsImageItem<T>::hoverMoveEvent(QGraphicsSceneHoverEvent *event) {
    if (_interaction) {
        _interaction->hoverMoved(this, event);
        return;
    }
    return QGraphicsPixmapItem::hoverMoveEvent(event);
}

template <class T>
void QGraphicsImageItem<T>::hoverLeaveEvent(QGraphicsSceneHoverEvent *event) {
    if (_interaction) {
        _interaction->hoverLeft(this, event);
        return;
    }
    return QGraphicsPixmapItem::hoverLeaveEvent(event);
}


#endif /* defined(__ParticleGuidedRegistration__QGraphicsImageItem__) */