//
//  QGraphicsEventItem.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/19/13.
//
//

#ifndef __ParticleGuidedRegistration__QGraphicsEventItem__
#define __ParticleGuidedRegistration__QGraphicsEventItem__

#include <iostream>
#include <QGraphicsItem>
#include <QGraphicsSceneMouseEvent>

template <class T>
class QGraphicsEventItem: public T {
public:
    class Listener {
        friend class QGraphicsEventItem;

    protected:
        virtual void mousePressed(QGraphicsEventItem<T>*, QGraphicsSceneMouseEvent*) = 0;
        virtual void mouseReleased(QGraphicsEventItem<T>*, QGraphicsSceneMouseEvent*) = 0;
    };

public:
    QGraphicsEventItem(): _listener(NULL) {}
    QGraphicsEventItem(QGraphicsItem* parent = NULL): T(parent), _listener(NULL) {
    }
    virtual ~QGraphicsEventItem() {}

    inline void setListener(Listener* listener) {
        _listener = listener;
    }

protected:
    void mousePressEvent(QGraphicsSceneMouseEvent* event);
    void mouseReleaseEvent(QGraphicsSceneMouseEvent* event);

private:
    Listener* _listener;
};

template <class T>
void QGraphicsEventItem<T>::mousePressEvent(QGraphicsSceneMouseEvent *event) {
    if (_listener) {
        _listener->mousePressed(this, event);
    }
}

template <class T>
void QGraphicsEventItem<T>::mouseReleaseEvent(QGraphicsSceneMouseEvent *event) {
    if (_listener) {
        _listener->mouseReleased(this, event);
    }
}


#endif /* defined(__ParticleGuidedRegistration__QGraphicsEventItem__) */
