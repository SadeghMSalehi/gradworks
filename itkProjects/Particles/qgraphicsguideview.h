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
#include "QGraphicsView"
#include <QList>


class QGraphicsGuideView: public QGraphicsView {
    Q_OBJECT
private:
    enum { DrawingMode, TransformMode } _currentMode;
    enum { Union, Subtract } _pathMergeMode;
    QPointF _drawingPressPoint;
    QPainterPath _drawingPath;
    QGraphicsPathItem* _drawingPathItem;
    QPainterPath _totalPath;
    QGraphicsPathItem* _totalPathItem;

    QList<QGraphicsPathItem*> _pathItemList;

public:
    QGraphicsGuideView(QWidget* parent = NULL);
    ~QGraphicsGuideView();

    void drawingMode();
    void transformMode();

    void mousePressEvent(QMouseEvent* event);
    void mouseMoveEvent(QMouseEvent* event);
    void mouseReleaseEvent(QMouseEvent* event);

    void drawingMousePressEvent(QMouseEvent* event);
    void drawingMouseMoveEvent(QMouseEvent* event);
    void drawingMouseReleaseEvent(QMouseEvent* event);


};

#endif /* defined(__ParticleGuidedRegistration__QGraphicsGuideView__) */
