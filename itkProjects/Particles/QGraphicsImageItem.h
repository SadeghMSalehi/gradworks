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
#include <qgraphicsitem.h>

#include "piImageHistogram.h"


extern unsigned int __grays[256];

template <class T>
class QGraphicsImageItem: public QGraphicsPixmapItem {
public:
    enum Flags { NoFlip, LeftRight, UpDown };

    QGraphicsImageItem(QGraphicsItem* parent = NULL): QGraphicsPixmapItem(parent) {
        _range[0] = 0;
        _range[1] = 255;
    }
    virtual ~QGraphicsImageItem() {}

    T getRange(int i);
    void setRange(T min, T max);
    void setFlip(Flags flags);
    void refresh();

    void setImage(T* inputBuffer, int w, int h);
    template <class S> void setImage(typename S::Pointer image, bool computeRange = false);
    static void convertToGrayscale(T*, double, double, QImage&);

private:
    T* _inputBuffer;
    T _range[2];
    QImage _grayImage;
};

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
    refresh();
}


template <class T> template <class S>
void QGraphicsImageItem<T>::setImage(typename S::Pointer image, bool computeRange) {
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
        _range[0] = _range[1] = *buff;
        for (int i = 1; i < nelems; i++, buff++) {
            _range[0] = std::min(*buff, _range[0]);
            _range[1] = std::max(*buff, _range[1]);
        }
    }
    setImage(image->GetBufferPointer(), w, h);
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

    const float pixelRange = _range[1] - _range[0];
    const float pixelMin = _range[0];
    convertToGrayscale(_inputBuffer, pixelMin, pixelRange, _grayImage);

    QPixmap pixmap = QPixmap::fromImage(_grayImage);
    setPixmap(pixmap);
}

template <class T>
void QGraphicsImageItem<T>::convertToGrayscale(T* inputBuffer,
                                               const double min, const double range, QImage& outputImage) {
    const int width = outputImage.width();
    for (int i = 0; i < outputImage.height(); i++) {
        uchar* outputBuffer = outputImage.scanLine(i);
        for (int j = 0; j < width; j++, outputBuffer++, inputBuffer++) {
            *outputBuffer = uchar(256u * (*inputBuffer-min)/(range));
        }
    }
}

#endif /* defined(__ParticleGuidedRegistration__QGraphicsImageItem__) */