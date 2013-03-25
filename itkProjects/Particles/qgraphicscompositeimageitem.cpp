//
//  qgraphicscompositeimageitem.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/24/13.
//
//
#include "qgraphicscompositeimageitem.h"

#include "piImageDef.h"
#include "piImageIO.h"
#include "QGraphicsSceneMouseEvent"

using namespace pi;

typedef pi::ImageDisplay<RealImage> ImageDisplayType;

bool QGraphicsCompositeImageItem::CheckCompositeBuffer(int id1) {
    if (_imageDisplays == NULL) {
        return false;
    }

    ImageIO<RealImage> io;
    ImageDisplayType* img1 = _imageDisplays->at(id1);

    // check foreground image
    RealImage::Pointer fImg = img1->GetResampled(_resampleIdx);
    if (fImg.IsNull()) {
        return false;
    }

    // check composite image
    if (_compositeImage.IsNull()
        || _compositeImage->GetBufferedRegion() != fImg->GetBufferedRegion()) {
        _compositeImage = io.CopyImage(fImg);
    }
    return _compositeImage.IsNotNull();
}

void QGraphicsCompositeImageItem::CompositeAlpha(int id1, int id2) {
    if (!CheckCompositeBuffer(id1)) {
        return;
    }

    _viewMin = SHRT_MIN;
    _viewMax = SHRT_MAX;

    // prepare input and output buffer
    ImageDisplayType* img1 = _imageDisplays->at(id1);
    RealImage::Pointer fImg = img1->GetResampled(_resampleIdx);
    const int nElems = fImg->GetPixelContainer()->Size();
    DataReal* cBuf = _compositeImage->GetBufferPointer();
    DataReal* fBuf = fImg->GetBufferPointer();

    if (id2 >= 0) {
        // if there is more than one image,
        ImageDisplayType* img2 = _imageDisplays->at(id2);
        RealImage::Pointer mImg = _imageDisplays->at(_movingId)->GetResampled(_resampleIdx);
        DataReal* mBuf = mImg->GetBufferPointer();
        for (int i = 0; i < nElems; i++) {
            const DataReal f = Clamp(fBuf[i], img1->histogram.rangeMin, img1->histogram.rangeMax);
            const DataReal m = Clamp(mBuf[i], img2->histogram.rangeMin, img2->histogram.rangeMax);
            cBuf[i] = _alpha*f+(1-_alpha)*m;
        }
    } else {
        // for the case of single image
        for (int i = 0; i < nElems; i++) {
            const DataReal f = Clamp(fBuf[i], img1->histogram.rangeMin, img1->histogram.rangeMax);
            cBuf[i] = _alpha * f;
        }
    }
}

void QGraphicsCompositeImageItem::CompositeCheckerBoard(int id1, int id2) {
    if (!CheckCompositeBuffer(id1)) {
        return;
    }

    if (id2 < 0 || id2 >= _imageDisplays->Count()) {
        return;
    }
    _viewMin = SHRT_MIN;
    _viewMax = SHRT_MAX;
    
    // prepare input and output buffer
    ImageDisplayType* img1 = _imageDisplays->at(id1);
    RealImage::Pointer fImg = img1->GetResampled(_resampleIdx);
    DataReal* fBuf = fImg->GetBufferPointer();

    ImageDisplayType* img2 = _imageDisplays->at(id2);
    RealImage::Pointer mImg = _imageDisplays->at(_movingId)->GetResampled(_resampleIdx);
    DataReal* mBuf = mImg->GetBufferPointer();

    DataReal* cBuf = _compositeImage->GetBufferPointer();

    // FIXME: need to select appropriate size, probably function at imageDisplays
    int w = _imageDisplays->GetGrid(_resampleIdx).Width();
    int h = _imageDisplays->GetGrid(_resampleIdx).Height();

    const int cbszW = w / _cbCols;
    const int cbszH = h / _cbRows;

    cout << cbszW << endl;
    cout << cbszH << endl;


    int k = 0;
    for (int j = 0; j < h; j++) {
        for (int i = 0; i < w; i++, k++) {
            const DataReal f = Clamp(fBuf[k], img1->histogram.rangeMin, img1->histogram.rangeMax);
            const DataReal m = Clamp(mBuf[k], img2->histogram.rangeMin, img2->histogram.rangeMax);
            if ((i/cbszW)%2 == (j/cbszH)%2) {
                cBuf[k] = f;
            } else {
                cBuf[k] = m;
            }
        }
    }
}

void QGraphicsCompositeImageItem::CompositeDifference(int id1, int id2) {
    if (!CheckCompositeBuffer(id1)) {
        return;
    }

    if (id2 < 0 || id2 >= _imageDisplays->Count()) {
        return;
    }

    // prepare input and output buffer
    ImageDisplayType* img1 = _imageDisplays->at(id1);
    RealImage::Pointer fImg = img1->GetResampled(_resampleIdx);
    DataReal* fBuf = fImg->GetBufferPointer();

    ImageDisplayType* img2 = _imageDisplays->at(id2);
    RealImage::Pointer mImg = _imageDisplays->at(_movingId)->GetResampled(_resampleIdx);
    DataReal* mBuf = mImg->GetBufferPointer();

    DataReal* cBuf = _compositeImage->GetBufferPointer();

    int nElems = fImg->GetPixelContainer()->Size();
    for (int i = 0; i < nElems; i++) {
        const DataReal f = Clamp(fBuf[i], img1->histogram.rangeMin, img1->histogram.rangeMax);
        const DataReal m = Clamp(mBuf[i], img2->histogram.rangeMin, img2->histogram.rangeMax);
        cBuf[i] = std::abs(m-f);
        if (i == 0) {
            _viewMin = cBuf[i];
            _viewMax = cBuf[i];
        } else {
            _viewMin = std::min(_viewMin, cBuf[i]);
            _viewMax = std::max(_viewMax, cBuf[i]);
        }
    }
}

void QGraphicsCompositeImageItem::Refresh(int id1, int id2) {
    if (_imageDisplays == NULL || _imageDisplays->Count() <= id1 || id1 < 0) {
        return;
    }

    m_id1 = id1;
    m_id2 = id2;

    if (_compositionMode == Alpha) {
        CompositeAlpha(id1, id2);
    } else if (_compositionMode == CheckerBoard) {
        CompositeCheckerBoard(id1, id2);
    } else if (_compositionMode == IntensityDifference) {
        CompositeDifference(id1, id2);
    }

    // if the composition image is not generated,
    if (_compositeImage.IsNull()) {
        return;
    }

    // convert to color image
    // assume that the composite image ranges from 0 to 65535 (in short range)
    ColorFilterType::Pointer colorFilter = ColorFilterType::New();
    colorFilter->SetInput(_compositeImage);
    colorFilter->UseManualScalingOn();
    colorFilter->SetMinimumValue(_viewMin);
    colorFilter->SetMaximumValue(_viewMax);
    colorFilter->Update();
    _rgbImage = colorFilter->GetOutput();

    // may have to correct to deal with various slice
    int w = _imageDisplays->GetGrid(_resampleIdx).Width();
    int h = _imageDisplays->GetGrid(_resampleIdx).Height();

    setPixmap(QPixmap::fromImage(QImage((unsigned char*) _rgbImage->GetBufferPointer(), w, h, QImage::Format_ARGB32)));
    update();
}

void QGraphicsCompositeImageItem::mousePressEvent(QGraphicsSceneMouseEvent *event) {
    _mousePressedButton = event->button();
    if (_mousePressedButton == Qt::LeftButton) {
        _mousePressedPoint = event->pos();
        if (_imageDisplays->IsValidId(m_id2)) {
            _imageDisplays->at(m_id2)->SetOrigin();
            Refresh(m_id1, m_id2);
        }
    }
}

void QGraphicsCompositeImageItem::mouseMoveEvent(QGraphicsSceneMouseEvent *event) {
    if (_mousePressedButton == Qt::LeftButton) {
        QPointF pos = event->pos();
        QPointF delta = pos - _mousePressedPoint;

        if (_imageDisplays->IsValidId(m_id2)) {
            _imageDisplays->at(m_id2)->SetOriginFromResampling(_resampleIdx, delta.x(), delta.y());
            Refresh(m_id1, m_id2);
        }
    }
}


void QGraphicsCompositeImageItem::mouseReleaseEvent(QGraphicsSceneMouseEvent *event) {
    if (_mousePressedButton == Qt::LeftButton) {
        if (_imageDisplays->IsValidId(m_id2)) {
            _imageDisplays->at(m_id2)->SetOrigin();
            emit originChanged();
        }
    }
    _mousePressedButton = Qt::NoButton;
}