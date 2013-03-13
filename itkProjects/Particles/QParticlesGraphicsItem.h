//
//  QParticlesGraphicsItem.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/12/13.
//
//

#ifndef __ParticleGuidedRegistration__QParticlesGraphicsItem__
#define __ParticleGuidedRegistration__QParticlesGraphicsItem__

#include <iostream>

#include "QGraphicsItem"
#include "QPen"
#include "QBrush"
#include "QRectF"

#include "piParticleCore.h"

class QParticlesGraphicsItem : public QGraphicsItem {
public:
    QParticlesGraphicsItem();
    virtual ~QParticlesGraphicsItem();

    void SetPen(QPen pen);
    void SetBrush(QBrush brush);
    void SetParticles(pi::Particle* particles, int n);

    QRectF boundingRect() const;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);

private:
    pi::Particle* m_Particles;
    int m_Counts;
    
    QRectF m_Bounds;
    QPen m_Pen;
    QBrush m_Brush;
};


#endif /* defined(__ParticleGuidedRegistration__QParticlesGraphicsItem__) */