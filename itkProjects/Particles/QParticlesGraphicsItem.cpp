//
//  QParticlesGraphicsItem.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/12/13.
//
//

#include "QParticlesGraphicsItem.h"


QParticlesGraphicsItem::QParticlesGraphicsItem() {

}

QParticlesGraphicsItem::~QParticlesGraphicsItem() {

}


void QParticlesGraphicsItem::SetPen(QPen pen) {
    m_Pen = pen;
}

void QParticlesGraphicsItem::SetBrush(QBrush brush) {
    m_Brush = brush;
}

void QParticlesGraphicsItem::SetParticles(pi::Particle* particles, int n) {
    m_Particles = particles;
    m_Counts = n;
}

QRectF QParticlesGraphicsItem::boundingRect() const {
    return m_Bounds;
}

void QParticlesGraphicsItem::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget) {

}