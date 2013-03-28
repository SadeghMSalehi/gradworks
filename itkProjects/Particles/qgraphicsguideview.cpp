//
//  QGraphicsGuideView.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/24/13.
//
//

#include "qgraphicsguideview.h"
#include "QMouseEvent"
#include "QGraphicsScene"
#include <QPainterPath>
#include <QGraphicsPathItem>

using namespace std;

QGraphicsGuideView::QGraphicsGuideView(QWidget* parent): QGraphicsView(parent) {
    _currentMode = TransformMode;
    _pathMergeMode = Union;
    _totalPathItem = NULL;
}

QGraphicsGuideView::~QGraphicsGuideView() {
}

void QGraphicsGuideView::drawingMode() {
    _currentMode = DrawingMode;
}

void QGraphicsGuideView::transformMode() {
    _currentMode = TransformMode;
}

void QGraphicsGuideView::mousePressEvent(QMouseEvent* event) {
    if (DrawingMode == _currentMode) {
        drawingMousePressEvent(event);
    } else {
        QGraphicsView::mousePressEvent(event);
    }
}

void QGraphicsGuideView::mouseMoveEvent(QMouseEvent* event) {
    if (DrawingMode == _currentMode) {
        drawingMouseMoveEvent(event);
    } else {
        QGraphicsView::mouseMoveEvent(event);
    }
}

void QGraphicsGuideView::mouseReleaseEvent(QMouseEvent* event) {
    if (DrawingMode == _currentMode) {
        drawingMouseReleaseEvent(event);
    } else {
        QGraphicsView::mouseReleaseEvent(event);
    }
}

void QGraphicsGuideView::drawingMousePressEvent(QMouseEvent *event) {
    QPointF pos = mapToScene(event->pos());
    _drawingPressPoint = pos;
    _drawingPath = QPainterPath(_drawingPressPoint);

    if (this->scene() != NULL) {
        _drawingPathItem = this->scene()->addPath(_drawingPath);
        _drawingPathItem->setOpacity(0.2);
    }
    if (event->buttons() & Qt::RightButton) {
        _drawingPathItem->setPen(QPen(Qt::red,3));
        _pathMergeMode = Subtract;
    } else {
        _drawingPathItem->setPen(QPen(Qt::green,3));
        _pathMergeMode = Union;
    }
}

void QGraphicsGuideView::drawingMouseMoveEvent(QMouseEvent *event) {
    QPointF pos = mapToScene(event->pos());
    _drawingPath.lineTo(pos);
    _drawingPathItem->setPath(_drawingPath);
    _drawingPathItem->update();
}

void QGraphicsGuideView::drawingMouseReleaseEvent(QMouseEvent *event) {
    if (_pathMergeMode == Subtract) {
        _totalPath = _totalPath.subtracted(_drawingPath);
    } else {
        _totalPath = _totalPath.united(_drawingPath);
    }
    if (_totalPathItem == NULL) {
        _totalPathItem = this->scene()->addPath(_totalPath);
        _totalPathItem->setBrush(QBrush(Qt::yellow,Qt::SolidPattern));
        _totalPathItem->setOpacity(0.2);
    }
    _totalPathItem->setPath(_totalPath);
    _totalPathItem->setPen(Qt::NoPen);
    _totalPathItem->update();
    this->scene()->removeItem(_drawingPathItem);
}