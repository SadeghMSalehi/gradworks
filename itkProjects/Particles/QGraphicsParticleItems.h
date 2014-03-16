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
#include "vnlCommon.h"

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
    void setReferenceSubject(pi::ParticleSubject* refSubject);
    
    void useScalars(bool);
    void setScalars(VNLVector scalars);


    void createParticles(pi::ParticleSubject* subject);
    void updateParticles();
    void clearParticles();
    void useParticleX(bool value);
    QGraphicsEllipseEventItem* getItem(int i);

    // be careful about this behavior
    void selectParticle(int);
    int getSelectedParticleId();

protected:
private:
    bool _useParticleX;
    bool _isParticleSelected;
    int _particleSelectedId;
    bool _isHide;
    QGraphicsScene* _scene;
    QGraphicsItem* _parentItem;
    pi::ParticleSubject* _subject;
    pi::ParticleSubject* _refSubject;
    QVector<QGraphicsEllipseEventItem*> _particleItems;
    QGraphicsEllipseEventItem::Listener* _listener;

    bool _useScalars;
    VNLVector _scalars;
};

#endif /* defined(__ParticleGuidedRegistration__QGraphicsParticleItems__) */