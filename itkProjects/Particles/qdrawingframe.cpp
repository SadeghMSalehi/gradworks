//
//  qdrawingframe.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/3/13.
//
//

#include "qdrawingframe.h"
#include <QPainter>
#include <QMouseEvent>

QDrawingFrame::QDrawingFrame(QWidget* widget, Qt::WindowType windowType): QFrame(widget, windowType) {

}

QDrawingFrame::~QDrawingFrame() {
    
}

void QDrawingFrame::paintEvent(QPaintEvent* event) {
    QPainter painter(this);
    painter.setPen(QPen(Qt::black));
    painter.setBrush(QBrush(Qt::black, Qt::SolidPattern));
    for (int i = 0; i < _drawingPath.elementCount(); i++) {
        QPainterPath::Element elem = _drawingPath.elementAt(i);
        painter.drawEllipse(elem.x-1, elem.y-1, 3, 3);
    }
    painter.setBrush(Qt::NoBrush);
    painter.drawPath(_drawingPath);
}

void QDrawingFrame::mousePressEvent(QMouseEvent* event) {
    _drawingPath.lineTo(event->pos());
}

void QDrawingFrame::mouseMoveEvent(QMouseEvent* event) {
    _drawingPath.lineTo(event->pos());
}

void QDrawingFrame::mouseReleaseEvent(QMouseEvent* event) {
    update();
}
