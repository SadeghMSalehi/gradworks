//
//  myQGraphicsView.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 11/19/12.
//
//

#include "myQGraphicsView.h"
#include "QMouseEvent"

void myQGraphicsView::mousePressEvent(QMouseEvent* event) {
    emit mousePressed(event);
}

void myQGraphicsView::mouseReleaseEvent(QMouseEvent* event) {
    emit mouseReleased(event);
}

void myQGraphicsView::mouseMoveEvent(QMouseEvent *event) {
    emit mouseMoved(event);
}