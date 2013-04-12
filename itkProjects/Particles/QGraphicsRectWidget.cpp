//
//  QGraphicsRectWidget.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/9/13.
//
//

#include "QGraphicsRectWidget.h"
#include <QPainter>
#include <QCursor>
#include <QGraphicsSceneHoverEvent>

using namespace std;

QGraphicsRectWidget::QGraphicsRectWidget(QGraphicsItem* parent): QAbstractGraphicsShapeItem(parent) {
    _hovering = false;
    _hoverPen = QPen(Qt::red, 1, Qt::SolidLine);
    _hoverPen.setCosmetic(true);
    _normalPen = QPen(Qt::green, 1, Qt::SolidLine);
    _normalPen.setCosmetic(true);

    setZValue(100);
    setAcceptHoverEvents(true);
    setAcceptedMouseButtons(Qt::LeftButton);

    _rect = QRectF(0, 0, 0, 0);
}

QGraphicsRectWidget::~QGraphicsRectWidget() {

}

void QGraphicsRectWidget::setSize(int n) {
    prepareGeometryChange();

    _rect.setWidth(n);
    _rect.setHeight(n);

    update();
}

QRectF QGraphicsRectWidget::boundingRect() const {
    return _rect;
}

void QGraphicsRectWidget::paint(QPainter *painter,
                                const QStyleOptionGraphicsItem *option, QWidget *widget) {
    QRectF rect = _rect;
    if (_hovering) {
        painter->setPen(_hoverPen);
    } else {
        painter->setPen(_normalPen);
    }
    painter->drawRect(rect);
}

void QGraphicsRectWidget::moveWidget(QPointF pos) {
    this->setPos(pos);
}

void QGraphicsRectWidget::hoverEnterEvent(QGraphicsSceneHoverEvent *e) {
    _hovering = true;
    update();
}

void QGraphicsRectWidget::hoverMoveEvent(QGraphicsSceneHoverEvent *e) {

}

void QGraphicsRectWidget::hoverLeaveEvent(QGraphicsSceneHoverEvent *e) {
    _hovering = false;
    update();
}

void QGraphicsRectWidget::mousePressEvent(QGraphicsSceneMouseEvent *e) {
    _startingPos = mapToItem(this, e->pos());
}

void QGraphicsRectWidget::mouseMoveEvent(QGraphicsSceneMouseEvent *e) {
    QPointF pos = mapToScene(e->pos());
    QPoint newPos(pos.x()-_startingPos.x()+.5, pos.y()-_startingPos.y()+.5);
    setPos(newPos);
}

void QGraphicsRectWidget::mouseReleaseEvent(QGraphicsSceneMouseEvent *e) {
    QPointF pos = mapToScene(e->pos());
    QPoint newPos(pos.x()-_startingPos.x()+.5, pos.y()-_startingPos.y()+.5);
    setPos(newPos);
    emit widgetMoved(pos);
//    QPointF pos = e->pos();
//    cout << pos.x() << "," << pos.y() << endl;
}