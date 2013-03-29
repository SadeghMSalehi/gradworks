//
//  QGraphicsGuideView.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/24/13.
//
//

#ifndef __ParticleGuidedRegistration__QGraphicsGuideView__
#define __ParticleGuidedRegistration__QGraphicsGuideView__

#include <iostream>
#include "qtypedef.h"
#include <QGraphicsView>
#include <QBitmap>
#include "piImageSlice.h"

class QAbstractGraphicsShapeItem;

class QGraphicsGuideView: public QGraphicsView {
    Q_OBJECT
public:
    enum ToolType { FreePath, FreeBrush, EraseBrush };
    
private:
    enum { DrawingMode, TransformMode } _currentMode;
    ToolType _currentTool;

    QPointF _drawingPressPoint;
    QPainterPath _drawingPath;
    QGraphicsPathItem* _drawingPathItem;

    QList<QGraphicsPathItem*> _pathItemList;

    // temporary canvas for image drawing
    uchar _foregroundLabel;
    uchar _backgroundLabel;
    QImage _drawingCanvas;

    // visualization of current label
    QSize _canvasSize;
    QImage _labelImage;
    QGraphicsPixmapItem* _labelItem;

    // brush cursor
    qreal _brushRadius;
    QAbstractGraphicsShapeItem* _brushCursor;

    int _volumeSliceIdx;
    pi::SliceDirectionEnum _volumeSliceDirection;
    pi::AIRLabel::Pointer _volumeSlice;
    pi::AIRLabel::Pointer _labelVolume;


private: // private methods
    void paintBrushEllipse(QRectF&);
    void canvasToSlice();
    void updateLabelItem();

    bool isLabelLoaded();
    void sliceToVolume();
    void volumeToSlice();
    void sliceToLabelImage();

    inline uint colorTable(pi::AIRClass p) {
        QColor c;
        static uint colors[5] = { 0, qRgb(0xff,0xff,0), qRgb(0,0xff,0), qRgb(0xff,0,0xff), qRgb(0,0xff,0xff) };
        return colors[p];
    }

public:
    QGraphicsGuideView(QWidget* parent = NULL);
    ~QGraphicsGuideView();

    void drawingMode(ToolType tool);
    void transformMode();

    void mousePressEvent(QMouseEvent* event);
    void mouseMoveEvent(QMouseEvent* event);
    void mouseReleaseEvent(QMouseEvent* event);

    void pathStartEvent(QMouseEvent* event);
    void pathMoveEvent(QMouseEvent* event);
    void pathEndEvent(QMouseEvent* event);

    void brushStartEvent(QMouseEvent* event);
    void brushMoveEvent(QMouseEvent* event);
    void brushEndEvent(QMouseEvent* event);

public slots:
    void segmentationCleared();
    void labelOpacityChanged(int n);
    void sliceChanged(pi::SliceDirectionEnum dir, int n);
    void labelVolumeChanged(pi::AIRLabel::Pointer);
    void createLabelVolume(pi::AIRImage::Pointer);
    void saveLabelVolume(QString& fileName);
    void propagateLabel(pi::IntVector& targetSlices);
};

#endif /* defined(__ParticleGuidedRegistration__QGraphicsGuideView__) */
