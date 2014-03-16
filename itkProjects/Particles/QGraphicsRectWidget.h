//
//  QGraphicsRectWidget.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/9/13.
//
//

#ifndef __ParticleGuidedRegistration__QGraphicsRectWidget__
#define __ParticleGuidedRegistration__QGraphicsRectWidget__

#include <iostream>
#include <QAbstractGraphicsShapeItem>
#include <QPen>

class QGraphicsRectWidget: public QObject, public QAbstractGraphicsShapeItem {
    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)
    
public:
    QGraphicsRectWidget(QGraphicsItem* parent = NULL);
    virtual ~QGraphicsRectWidget();

    void setSize(int n);
    QRectF boundingRect() const;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget = NULL);
    void hoverEnterEvent(QGraphicsSceneHoverEvent* e);
    void hoverMoveEvent(QGraphicsSceneHoverEvent* e);
    void hoverLeaveEvent(QGraphicsSceneHoverEvent* e);
    void mousePressEvent(QGraphicsSceneMouseEvent* e);
    void mouseMoveEvent(QGraphicsSceneMouseEvent* e);
    void mouseReleaseEvent(QGraphicsSceneMouseEvent* e);

public slots:
    void moveWidget(QPointF pos);

signals:
    void widgetMoved(QPointF pos);

protected:
    QPointF _startingPos;
    QRectF _rect;
    bool _hovering;
    QPen _hoverPen;
    QPen _normalPen;
};
#endif /* defined(__ParticleGuidedRegistration__QGraphicsRectWidget__) */
