//
//  qmygraphicsview.cpp
//  imageParticles
//
//  Created by Joohwi Lee on 11/2/12.
//
//

#include "qmygraphicsview.h"
#include "QMouseEvent"
#include "iostream"


void QMyGraphicsView::mousePressEvent(QMouseEvent* event) {
    QPointF xy = this->mapToScene(event->pos());
    emit mousePressed(event);
}

void QMyGraphicsView::mouseReleaseEvent(QMouseEvent* event) {
    QPointF xy = this->mapToScene(event->pos());
    emit mouseReleased(event);
}