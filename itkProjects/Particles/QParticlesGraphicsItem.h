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
    enum ColorMode {
        Brush,
        Density,
        Label
    } colorMode;

    QParticlesGraphicsItem();
    virtual ~QParticlesGraphicsItem();

    void SetPen(QPen pen);
    void SetBrush(QBrush brush);
    void SetParticles(pi::Particle* particles, int n);
    
    void ColorModeToDensity() { colorMode = Density; }
    void ColorModeToLabel() { colorMode = Label; }
    void ColorModeToBrush() { colorMode = Brush; }

    QRectF boundingRect() const;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);

private:
    pi::Particle* m_Particles;
    int m_Counts;

    float _md[2];
    QRectF m_Bounds;
    QPen m_Pen;
    QBrush m_Brush;
};


#endif /* defined(__ParticleGuidedRegistration__QParticlesGraphicsItem__) */