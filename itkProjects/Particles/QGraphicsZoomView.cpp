//
//  QGraphicsZoomView.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 5/4/13.
//
//

#include "QGraphicsZoomView.h"
#include <QMouseEvent>
#include <QGraphicsScene>

using namespace std;

QGraphicsZoomView::QGraphicsZoomView(QWidget* parent): QGraphicsView(parent) {
    
}

QGraphicsZoomView::~QGraphicsZoomView() {
    
}

void QGraphicsZoomView::mousePressEvent(QMouseEvent *event) {
    if (event->buttons() & Qt::RightButton) {
        _trackingStartPoint = event->pos();
        QTransform xform = transform();
        _scale = xform.m11();
    }
}

void QGraphicsZoomView::mouseMoveEvent(QMouseEvent *event) {
    QPointF pos = event->pos();
    QPointF delta = pos - _trackingStartPoint;

    int scaleFactor = 50 - delta.y();
    scaleFactor = scaleFactor < 0.01 ? 0.01 : scaleFactor;
    float scale = _scale * (scaleFactor / 50.0);

    QTransform transform;
    transform.scale(scale, scale);

    setTransform(transform);
}

void QGraphicsZoomView::mouseReleaseEvent(QMouseEvent *event) {

}

void QGraphicsZoomView::enterEvent(QEvent *event) {

}

void QGraphicsZoomView::leaveEvent(QEvent *event) {

}