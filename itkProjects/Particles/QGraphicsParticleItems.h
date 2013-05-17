//
//  QGraphicsParticleItems.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/28/13.
//
//

#ifndef __ParticleGuidedRegistration__QGraphicsParticleItems__
#define __ParticleGuidedRegistration__QGraphicsParticleItems__

#include <iostream>
#include <QGraphicsEllipseItem>
#include "QGraphicsEventItem.h"

namespace pi {
    class ParticleSubject;
}

class QGraphicsScene;
typedef QGraphicsEventItem<QGraphicsEllipseItem> QGraphicsEllipseEventItem;

class QGraphicsParticleItems {
public:
    QGraphicsParticleItems();
    virtual ~QGraphicsParticleItems();

    void hideParticles(bool);
    void setListener(QGraphicsEllipseEventItem::Listener* listener);
    void setScene(QGraphicsScene* scene);
    void setParentItem(QGraphicsItem* parentItem);
    void createParticles(pi::ParticleSubject* subject);
    void updateParticles();
    void clearParticles();

    // be careful about this behavior
    void selectParticle(int);
    int getSelectedParticleId();

protected:
private:
    bool _isParticleSelected;
    int _particleSelectedId;
    bool _isHide;
    QGraphicsScene* _scene;
    QGraphicsItem* _parentItem;
    pi::ParticleSubject* _subject;
    QVector<QGraphicsEllipseEventItem*> _particleItems;
    QGraphicsEllipseEventItem::Listener* _listener;
};

#endif /* defined(__ParticleGuidedRegistration__QGraphicsParticleItems__) */