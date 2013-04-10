//
//  myQGraphicsView.cpp
//  ParticlesGUI
//
//  Created by Joohwi Lee on 11/19/12.
//
//

#include "myQGraphicsView.h"
#include "QMouseEvent"

void myQGraphicsView::enterEvent(QEvent *event) {

}

void myQGraphicsView::leaveEvent(QEvent *event) {
    _rightMode = None;
}

void myQGraphicsView::mousePressEvent(QMouseEvent* event) {
    if (event->button() == Qt::RightButton) {
        _rightMode = Zooming;
        _buttonPressPos = event->pos();
        _buttonPressViewRect = contentsRect();
    }
    QGraphicsView::mousePressEvent(event);
}

void myQGraphicsView::mouseReleaseEvent(QMouseEvent* event) {
    if (event->button() == Qt::RightButton) {
        _rightMode = None;
    }
    QGraphicsView::mouseReleaseEvent(event);
}


void myQGraphicsView::mouseMoveEvent(QMouseEvent *event) {
    if (_rightMode == Zooming) {
        QPoint drag = event->pos() - _buttonPressPos;
        QRect rect = _buttonPressViewRect;
//        std::cout << drag.y() << std::endl;
        rect.adjust(-drag.y(), -drag.y(), drag.y(), drag.y());
        QRectF viewRect = mapToScene(rect).boundingRect();
//        fitInView(viewRect, Qt::KeepAspectRatio);
    }

    QGraphicsView::mouseMoveEvent(event);
}