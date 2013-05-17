//
//  QGraphicsDirectionItem.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 5/13/13.
//
//

#ifndef __ParticleGuidedRegistration__QGraphicsDirectionItem__
#define __ParticleGuidedRegistration__QGraphicsDirectionItem__

#include <iostream>
#include <cmath>

#include <QGraphicsLineItem>

class QGraphicsDirectionItem: public QGraphicsLineItem {
public:
    QGraphicsDirectionItem(QGraphicsItem* parentItem): QGraphicsLineItem(parentItem) {}
    virtual ~QGraphicsDirectionItem() {}

    void setDirection(QPointF direction) {
        _direction = direction;
        if (parentItem() != NULL) {
            QRectF rect = parentItem()->boundingRect();
            QPointF startingPos(rect.width() / 2, rect.height() / 2);

            int length = std::min(rect.width() / 2.0, rect.height() / 2.0);
            float mag = std::sqrt(direction.x()*direction.x() + direction.y()*direction.y());
            _direction *= (length/mag);

            QPointF endingPos = startingPos + _direction;
            this->setLine(QLineF(startingPos, endingPos));
        }
    }

private:
    QPointF _direction;
};
#endif /* defined(__ParticleGuidedRegistration__QGraphicsDirectionItem__) */
