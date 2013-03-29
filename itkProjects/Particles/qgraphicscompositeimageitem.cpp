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
#include "QGraphicsScene"

using namespace pi;

typedef pi::ImageDisplay<AIRImage> ImageDisplayType;

bool QGraphicsCompositeImageItem::CheckCompositeBuffer(int id1) {
    if (_imageDisplays == NULL) {
        return false;
    }

    ImageIO<AIRImage> io;
    ImageDisplayType* img1 = _imageDisplays->at(id1);

    // check foreground image
    AIRImage::Pointer fImg = img1->GetResampled(_resampleIdx);
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

    _viewMin = 0;
    _viewMax = USHRT_MAX;

    // prepare input and output buffer
    ImageDisplayType* img1 = _imageDisplays->at(id1);
    AIRImage::Pointer fImg = img1->GetResampled(_resampleIdx);
    const int nElems = fImg->GetPixelContainer()->Size();
    AIRPixel* cBuf = _compositeImage->GetBufferPointer();
    AIRPixel* fBuf = fImg->GetBufferPointer();
    UShortClamper fClamp(img1->histogram.rangeMin, img1->histogram.rangeMax);

    if (_imageDisplays->IsValidId(id2)) {
        // if there is more than one image,
        ImageDisplayType* img2 = _imageDisplays->at(id2);
        AIRImage::Pointer mImg = img2->GetResampled(_resampleIdx);
        AIRPixel* mBuf = mImg->GetBufferPointer();

        UShortClamper mClamp(img2->histogram.rangeMin, img2->histogram.rangeMax);
        for (int i = 0; i < nElems; i++) {
            const AIRPixel f = fClamp(fBuf[i]);
            const AIRPixel m = mClamp(mBuf[i]);
            cBuf[i] = AIRPixel(_alpha*f+(1-_alpha)*m);
        }
    } else {
        // for the case of single image
        for (int i = 0; i < nElems; i++) {
            const AIRPixel f = fClamp(fBuf[i]);
            cBuf[i] = AIRPixel(_alpha * f);
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
    _viewMin = 0;
    _viewMax = USHRT_MAX;
    
    // prepare input and output buffer
    ImageDisplayType* img1 = _imageDisplays->at(id1);
    AIRImage::Pointer fImg = img1->GetResampled(_resampleIdx);
    AIRPixel* fBuf = fImg->GetBufferPointer();

    ImageDisplayType* img2 = _imageDisplays->at(id2);
    AIRImage::Pointer mImg = img2->GetResampled(_resampleIdx);
    AIRPixel* mBuf = mImg->GetBufferPointer();

    AIRPixel* cBuf = _compositeImage->GetBufferPointer();

    // FIXME: need to select appropriate size, probably function at imageDisplays
    int w = _imageDisplays->GetGrid(_resampleIdx).Width();
    int h = _imageDisplays->GetGrid(_resampleIdx).Height();

    const int cbszW = w / _cbCols;
    const int cbszH = h / _cbRows;

    cout << cbszW << endl;
    cout << cbszH << endl;

    UShortClamper fClamp(img1->histogram.rangeMin, img1->histogram.rangeMax);
    UShortClamper mClamp(img2->histogram.rangeMin, img2->histogram.rangeMax);

    int k = 0;
    for (int j = 0; j < h; j++) {
        for (int i = 0; i < w; i++, k++) {
            const AIRPixel f = fClamp(fBuf[k]);
            const AIRPixel m = mClamp(mBuf[k]);
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
    AIRImage::Pointer fImg = img1->GetResampled(_resampleIdx);
    AIRPixel* fBuf = fImg->GetBufferPointer();

    ImageDisplayType* img2 = _imageDisplays->at(id2);
    AIRImage::Pointer mImg = img2->GetResampled(_resampleIdx);
    AIRPixel* mBuf = mImg->GetBufferPointer();
    AIRPixel* cBuf = _compositeImage->GetBufferPointer();


    UShortClamper fClamp(img1->histogram.rangeMin, img1->histogram.rangeMax);
    UShortClamper mClamp(img2->histogram.rangeMin, img2->histogram.rangeMax);

    int nElems = fImg->GetPixelContainer()->Size();
    for (int i = 0; i < nElems; i++) {
        const AIRPixel f = fClamp(fBuf[i]);
        const AIRPixel m = mClamp(mBuf[i]);
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
    typedef itk::ScalarToARGBColormapImageFilter<AIRImage, RGBAVolumeType> ColorFilterType;

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

    QSize size(w,h);
    if (m_imageSize != size) {
        m_imageSize = size;
        emit imageSizeChanged(m_imageSize);
    }
    update();
}

void QGraphicsCompositeImageItem::SetInteractionModeToNone() {
    QGraphicsScene* scene = this->scene();
    scene->removeItem(m_drawingImageItem);
    m_drawingImageItem = NULL;
}

void QGraphicsCompositeImageItem::SetInteractionModeToDrawing() {
    _mode = Drawing;

    QGraphicsScene* scene = this->scene();
    QRectF rect = this->boundingRect();

    if (rect.width() > 0 && rect.height() > 0) {
        m_drawingImage = QImage(rect.width(), rect.height(), QImage::Format_RGB16);
        m_drawingImageItem = scene->addPixmap(QPixmap::fromImage(m_drawingImage));
    }
}


void QGraphicsCompositeImageItem::mousePressEvent(QGraphicsSceneMouseEvent *event) {
    _mousePressedButton = event->button();
    if (_mode == None) {
        if (_mousePressedButton == Qt::LeftButton) {
            _mousePressedPoint = event->pos();
            if (_imageDisplays->IsValidId(m_id2)) {
                _imageDisplays->at(m_id2)->SetAffineTranslation();
                Refresh(m_id1, m_id2);
            }
        }
    } else if (_mode == Drawing) {
        drawingBeginEvent(event);
    }
}

void QGraphicsCompositeImageItem::mouseMoveEvent(QGraphicsSceneMouseEvent *event) {
    if (_mousePressedButton == Qt::LeftButton) {
        if (_mode == None) {
            QPointF pos = event->pos();
            QPointF delta = _mousePressedPoint - pos;

            if (_imageDisplays->IsValidId(m_id2)) {
                _imageDisplays->at(m_id2)->SetAffineTranslation(_resampleIdx, delta.x(), delta.y());
                emit translationChanged();
                Refresh(m_id1, m_id2);
            }
        } else if (_mode == Drawing) {
            drawingMoveEvent(event);
        }
    }
}


void QGraphicsCompositeImageItem::mouseReleaseEvent(QGraphicsSceneMouseEvent *event) {
    if (_mousePressedButton == Qt::LeftButton) {
        if (_mode == None) {
            if (_imageDisplays->IsValidId(m_id2)) {
                _imageDisplays->at(m_id2)->SetAffineTranslation();
            }
        } else {
            drawingFinishEvent(event);
        }
    }
    _mousePressedButton = Qt::NoButton;
}

void QGraphicsCompositeImageItem::drawingBeginEvent(QGraphicsSceneMouseEvent *event) {

}

void QGraphicsCompositeImageItem::drawingMoveEvent(QGraphicsSceneMouseEvent *event) {

}

void QGraphicsCompositeImageItem::drawingFinishEvent(QGraphicsSceneMouseEvent *event) {

}
