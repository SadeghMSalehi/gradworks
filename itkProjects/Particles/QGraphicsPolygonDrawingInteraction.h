//
//  QGraphicsPolygonDrawingInteraction.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 5/3/13.
//
//

#ifndef __ParticleGuidedRegistration__QGraphicsPolygonDrawingInteraction__
#define __ParticleGuidedRegistration__QGraphicsPolygonDrawingInteraction__

#include <iostream>
#include <QVector>
#include <QGraphicsPolygonItem>
#include <QPolygonF>
#include "QGraphicsImageItem.h"


#include "piFitCurve.h"


template <class T>
class QGraphicsPolygonDrawingInteraction: public QGraphicsItemInteraction<T> {
private:
    bool _interactionBegin;
    QPointF _startingPos;
    QPolygonF _points;
    QGraphicsPolygonItem* _polygonItem;
    pi::ParticleVector _particles;

public:
    QGraphicsPolygonDrawingInteraction() {
        static QPen pen(Qt::yellow, 1);
        pen.setCosmetic(false);

        _interactionBegin = false;
        _polygonItem = new QGraphicsPolygonItem();
        _polygonItem->setPen(pen);
        _polygonItem->setBrush(Qt::NoBrush);
        _polygonItem->setOpacity(0.5);
    }

    void reset() {
        _interactionBegin = false;
        _polygonItem->setParentItem(NULL);
        _polygonItem->hide();
        _points.clear();
    }

    QPolygonF getPoints() {
        return _points;
    }

    pi::ParticleVector& getParticles() {
        return _particles;
    }


    void setParticles(T* item, pi::ParticleVector& particles) {
        _particles = particles;
        _points.clear();

        for (int i = 0; i < _particles.size(); i++) {
            _points.append(QPointF(_particles[i].x[0], _particles[i].x[1]));
        }
        _polygonItem->setPolygon(_points);
        _polygonItem->setParentItem(item);
        _polygonItem->show();
    }


    virtual void mousePressed(T *item, QGraphicsSceneMouseEvent *event) {
        if (event->buttons() & Qt::LeftButton) {
            if (_points.size() == 0) {
                _polygonItem->setParentItem(item);
                _polygonItem->setPolygon(_points);
                _polygonItem->show();
            }
            _interactionBegin = true;
        } else {
            _points.clear();
            _polygonItem->setParentItem(NULL);
            _polygonItem->hide();
        }
    }

    virtual void mouseMoved(T *item, QGraphicsSceneMouseEvent *event) {
        if (_interactionBegin) {
            QPointF pos = event->pos();
            std::cout << pos.x() << ", " << pos.y() << std::endl;
            _points.append(event->pos());
            _polygonItem->setPolygon(_points);
        }
    }

    virtual void mouseReleased(T *item, QGraphicsSceneMouseEvent *event) {
        _interactionBegin = false;

        using namespace pi;

        CurveFitting curveFitting;
        ParticleVector particles;
        particles.resize(_points.size());

        for (int i = 0; i < particles.size(); i++) {
            particles[i].x[0] = _points[i].x();
            particles[i].x[1] = _points[i].y();
        }

        curveFitting.FitCurve(particles, true);
        _particles = curveFitting.GetResult();

        _points.clear();
        for (int i = 0; i < _particles.size(); i++) {
            _points.append(QPointF(_particles[i].x[0], _particles[i].x[1]));
        }

        _polygonItem->setPolygon(_points);
    }

    virtual void hoverEntered(T *item, QGraphicsSceneHoverEvent *event) {

    }

    virtual void hoverMoved(T *item, QGraphicsSceneHoverEvent *event) {

    }

    virtual void hoverLeft(T *item, QGraphicsSceneHoverEvent *event) {

    }

};

#endif /* defined(__ParticleGuidedRegistration__QGraphicsPolygonDrawingInteraction__) */