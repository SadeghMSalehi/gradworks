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
#include <itkRGBAPixel.h>
#include <itkImage.h>
#include <qgraphicsitem.h>

#include "piImageHistogram.h"
#include "itkScalarToARGBColormapImageFilter.h"

template <class T>
class QGraphicsImageItem: public QGraphicsPixmapItem {
public:
    enum { ImageDimension = T::ImageDimension };
    enum Flags { NoFlip, LeftRight, UpDown };
    typedef typename T::Pointer ImagePointer;
    typedef typename T::PixelType ImagePixel;
    typedef itk::RGBAPixel<unsigned char> ColorPixel;
    typedef itk::Image<ColorPixel, ImageDimension> ColorImage;
    typedef typename ColorImage::Pointer ColorImagePointer;
    typedef itk::ScalarToARGBColormapImageFilter<T, ColorImage> ColorFilter;
    typedef pi::ImageHistogram<T> HistogramType;

    QGraphicsImageItem(QGraphicsItem* parent = NULL): QGraphicsPixmapItem(parent) {
    }
    virtual ~QGraphicsImageItem() {}

    void setImage(ImagePointer image);
    void setImage(ImagePointer image, ImagePixel min, ImagePixel max);
    void setIntensityRange(ImagePixel min, ImagePixel max);
    void setFlip(Flags flags);
    void refresh();
    HistogramType& histogram();

private:
    ImagePointer _image;
    ColorImagePointer _colorImage;
    HistogramType _histogram;
    int _width;
    int _height;

    void computeSize();
};

template <class T>
typename QGraphicsImageItem<T>::HistogramType& QGraphicsImageItem<T>::histogram() {
    return _histogram;
}

template <class T>
void QGraphicsImageItem<T>::setImage(ImagePointer image) {
    if (image.IsNull()) {
        _image = image;
        _width = 0;
        _height = 0;
        setPixmap(QPixmap());
        return;
    }
    _image = image;
    _histogram.SetImage(image);
    computeSize();
}

template <class T>
void QGraphicsImageItem<T>::setImage(ImagePointer image, ImagePixel min, ImagePixel max) {
    if (image.IsNull()) {
        _image = image;
        _width = 0;
        _height = 0;
        setPixmap(QPixmap());
        return;
    }
    _image = image;
    _histogram.SetImage(image, false);
    _histogram.dataMin = _histogram.rangeMin = min;
    _histogram.dataMax = _histogram.rangeMax = max;
    computeSize();
}


template <class T>
void QGraphicsImageItem<T>::setIntensityRange(ImagePixel min, ImagePixel max) {
    _histogram.rangeMin = min;
    _histogram.rangeMax = max;
}

template <class T>
void QGraphicsImageItem<T>::setFlip(Flags flags) {
    if (flags == LeftRight) {
        QTransform transform;
        transform.setMatrix(-1, 0, 0, 0, 1, 0, _width, 0, 1);
        this->setTransform(transform);
    } else if (flags == UpDown) {
        QTransform transform;
        transform.setMatrix(1, 0, 0, 0, -1, 0, 0, _height, 1);
        this->setTransform(transform);
    } else {
        QTransform transform;
        transform.reset();
        this->setTransform(transform);
    }
}

template <class T>
void QGraphicsImageItem<T>::refresh() {
    typename ColorFilter::Pointer filter = ColorFilter::New();
    filter->SetInput(_image);
    filter->UseManualScalingOn();
    filter->SetMinimumValue(_histogram.dataMin);
    filter->SetMaximumValue(_histogram.dataMax);
    filter->Update();
    _colorImage = filter->GetOutput();
    _colorImage->DisconnectPipeline();

    QPixmap pixmap = QPixmap::fromImage(QImage((unsigned char*) _colorImage->GetBufferPointer(), _width, _height, QImage::Format_ARGB32));
    setPixmap(pixmap);
}

template <class T>
void QGraphicsImageItem<T>::computeSize() {
    typename T::RegionType region = _image->GetBufferedRegion();
    if (region.GetSize(2) < 2 || ImageDimension == 2) {
        _width = region.GetSize(0);
        _height = region.GetSize(1);
    } else if (region.GetSize(1) < 2) {
        _width = region.GetSize(0);
        _height = region.GetSize(2);
    } else if (region.GetSize(0) < 2) {
        _width = region.GetSize(1);
        _height = region.GetSize(2);
    }
}
#endif /* defined(__ParticleGuidedRegistration__QGraphicsImageItem__) */