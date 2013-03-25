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

class QGraphicsCompositeImageItem: public QObject, public QGraphicsPixmapItem {
    Q_OBJECT

private:
    enum CompositionMode { Alpha, CheckerBoard, IntensityDifference };

signals:
    void originChanged();

private:
    int _resampleIdx;
    double _alpha;
    int _cbRows, _cbCols;
    CompositionMode _compositionMode;
    pi::DataReal _viewMin, _viewMax;

    // fixed and moving image
    int m_id1, m_id2;
    QPointF _mousePressedPoint;
    Qt::MouseButton _mousePressedButton;

    // memory holder for qpixmap
    pi::RGBAVolumeType::Pointer _rgbImage;

    // memory holder for composite image
    pi::RealImage::Pointer _compositeImage;


    
private:
    typedef pi::ImageDisplayCollection<pi::RealImage> ImageCollectionType;
    ImageCollectionType* _imageDisplays;

    int _fixedId;
    int _movingId;

    // image displays
    bool CheckCompositeBuffer(int id1);
    void CompositeAlpha(int id1, int id2);
    void CompositeCheckerBoard(int id1, int id2);
    void CompositeDifference(int id1, int id2);

    inline pi::DataReal Clamp(pi::DataReal f, const pi::DataReal fMin, const pi::DataReal fMax) {
        if (f < fMin) {
            f = fMin;
        } else if (f > fMax) {
            f = fMax;
        }
        f = (f-fMin)/(fMax-fMin)*(SHRT_MAX)+SHRT_MIN;
        return f;
    }

public:
    QGraphicsCompositeImageItem(QGraphicsItem* parent = NULL): QGraphicsPixmapItem(parent), _resampleIdx(0) {
        _fixedId = 0;
        _movingId = 1;
        _imageDisplays = NULL;
        _compositionMode = Alpha;
        _alpha = 1;
        _cbRows = 4;
        _cbCols = 4;
    }

    void SetImageDisplays(typename pi::ImageDisplayCollection<pi::RealImage>* displays) {
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
};

#endif /* defined(__ParticleGuidedRegistration__qgraphicscompositeimageitem__) */
