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

extern uint __maxLabelColors;
extern uint __labelColors[254];
extern uint* __labelColorsPtr;

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

    // visualization of current label
    QRectF _updateRect;
    QSize _canvasSize;
    QImage _userDrawingCanvas;
    QImage _sliceImage;
    QImage _labelImage;
    QGraphicsPixmapItem* _labelItem;

    // drawing labels
    uchar _foregroundLabel;
    uchar _backgroundLabel;
    
    // brush cursor
    qreal _brushRadius;
    QAbstractGraphicsShapeItem* _brushCursor;
    
    int _volumeSliceIdx;
    pi::SliceDirectionEnum _volumeSliceDirection;
    pi::AIRLabelSlice _volumeSlice;
    pi::AIRLabel::Pointer _labelVolume;
    pi::AIRLabel::Pointer _copiedSlice;


private: // private methods
    void paintBrushEllipse(QRectF&);
    void userDrawingToSlice();
    void updateLabelItem();

    bool isLabelLoaded();
    void sliceToVolume(bool);
    void volumeToSlice();
    void sliceToLabelImage(bool);


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

    void setInteractionLabels(int f, int b);

    inline uint colorTable(pi::AIRClass p) {
        QColor c;
        return __labelColorsPtr[p];
    }

    pi::AIRLabelSlice getLabelSlice();
    pi::AIRLabel::Pointer getLabelVolume();

public slots:
    void segmentationCleared();
    void labelOpacityChanged(int n);
    void sliceChanged(pi::SliceDirectionEnum dir, int n);
    void labelVolumeChanged(pi::AIRLabel::Pointer);
    void createLabelVolumeIfNecessary(pi::AIRImage::Pointer);
    void loadLabelVolume(QString& fileName, pi::AIRImage::Pointer);
    void saveLabelVolume(QString& fileName);
    void propagateLabel(pi::IntVector& targetSlices);

    void copyLabel();
    void pasteLabel();
};

#endif /* defined(__ParticleGuidedRegistration__QGraphicsGuideView__) */
