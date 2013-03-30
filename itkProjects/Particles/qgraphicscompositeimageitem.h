//
//  qgraphicscompositeimageitem.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/24/13.
//
//

#ifndef __ParticleGuidedRegistration__qgraphicscompositeimageitem__
#define __ParticleGuidedRegistration__qgraphicscompositeimageitem__

#include <iostream>
#include <climits>
#include "piImageSlice.h"
#include "QGraphicsPixmapItem"

class UShortClamper {
public:
    const double Min;
    const double Max;

    UShortClamper(double min, double max): Min(min), Max(max) {
    }
    
    inline unsigned short operator()(double f) {
        return ushort((f-Min)/(Max-Min)*(USHRT_MAX));
    }
};

class FloatClamper {
public:
    const float Min;
    const float Max;

    FloatClamper(float min, float max): Min(min), Max(max) {
    }

    inline float operator()(float f) {
        return (f-Min)/(Max-Min);
    }
};

typedef FloatClamper Clamper;

class QGraphicsCompositeImageItem: public QObject, public QGraphicsPixmapItem {
    Q_OBJECT

private:
    enum CompositionMode { Alpha, CheckerBoard, IntensityDifference };
    enum InteractionMode { None, Drawing };

signals:
    void translationChanged();
    void imageSizeChanged(QSize size);

private:
    // fixed and moving image
    int m_id1, m_id2;
    QSize m_imageSize;

    int _resampleIdx;
    double _alpha;
    int _cbRows, _cbCols;
    pi::DataReal _viewMin, _viewMax;

    InteractionMode _mode;
    CompositionMode _compositionMode;

    QPointF _mousePressedPoint;
    Qt::MouseButton _mousePressedButton;

    // memory holder for qpixmap
    pi::RGBAVolumeType::Pointer _rgbImage;

    // memory holder for composite image
    pi::RealImage::Pointer _compositeImage;

    // drawing buffer for current image
    QImage m_drawingImage;
    QGraphicsPixmapItem* m_drawingImageItem;

private:
    typedef pi::ImageDisplayCollection<pi::AIRImage> ImageCollectionType;
    ImageCollectionType* _imageDisplays;

    // image displays
    bool CheckCompositeBuffer(int id1);
    void CompositeAlpha(int id1, int id2);
    void CompositeCheckerBoard(int id1, int id2);
    void CompositeDifference(int id1, int id2);


public:
    QGraphicsCompositeImageItem(QGraphicsItem* parent = NULL): QGraphicsPixmapItem(parent), _resampleIdx(0) {
        m_id1 = 0;
        m_id2 = 1;
        _imageDisplays = NULL;
        _compositionMode = Alpha;
        _alpha = 1;
        _cbRows = 4;
        _cbCols = 4;
        _mode = None;
        m_drawingImageItem = NULL;
    }

    QSize GetImageSize() {
        return m_imageSize;
    }

    void SetInteractionModeToNone();
    void SetInteractionModeToDrawing();

    void SetImageDisplays(ImageCollectionType* displays) {
        _imageDisplays = displays;
    }

    void CompositionModeToAlpha(double alpha) {
        _compositionMode = Alpha;
        _alpha = alpha;
    }

    void CompositionModeToCheckerBoard(int r, int c) {
        _compositionMode = CheckerBoard;
        _cbRows = r;
        _cbCols = c;
    }

    void CompositionMode(int option) {
        if (option == 0) {
            _compositionMode = IntensityDifference;
        }
    }

    void Refresh(int id1, int id2);


    void mousePressEvent(QGraphicsSceneMouseEvent* event);
    void mouseMoveEvent(QGraphicsSceneMouseEvent* event);
    void mouseReleaseEvent(QGraphicsSceneMouseEvent *event);

    void drawingBeginEvent(QGraphicsSceneMouseEvent* event);
    void drawingMoveEvent(QGraphicsSceneMouseEvent* event);
    void drawingFinishEvent(QGraphicsSceneMouseEvent *event);
};

#endif /* defined(__ParticleGuidedRegistration__qgraphicscompositeimageitem__) */
