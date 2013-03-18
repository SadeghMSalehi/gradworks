//
//  QParticlesGraphicsItem.cpp
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 3/12/13.
//
//

#include "QParticlesGraphicsItem.h"
#include "QPainter"
#include "QStyleOptionGraphicsItem"
#include "itkARGBColorFunction.h"
#include "piParticleCore.h"

QParticlesGraphicsItem::QParticlesGraphicsItem() {
    m_Counts = 0;
    m_Particles = NULL;
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
    float mx[2] = { 0 }, my[2] = { 0 };
    for (int i = 0; i < n; i++) {
        if (i == 0) {
            mx[0] = mx[1] = m_Particles[i].x[0];
            my[0] = my[1] = m_Particles[i].x[1];
            _md[0] = _md[1] = m_Particles[i].density;
        } else {
            mx[0] = ::min(m_Particles[i].x[0], mx[0]);
            mx[1] = ::max(m_Particles[i].x[0], mx[1]);
            my[0] = ::min(m_Particles[i].x[1], my[0]);
            my[1] = ::max(m_Particles[i].x[1], my[1]);
            _md[0] = ::min(m_Particles[i].density, _md[0]);
            _md[1] = ::max(m_Particles[i].density, _md[1]);
        }
    }
    m_Bounds = QRectF(mx[0]-3, my[0]-3, mx[1]-mx[0]+6, my[1]-my[0]+6);
}

QRectF QParticlesGraphicsItem::boundingRect() const {
    return m_Bounds;
}

void QParticlesGraphicsItem::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget) {
    if (m_Particles == NULL) {
        return;
    }
    typedef itk::RGBAPixel<unsigned char> RGBA;
    typedef itk::Function::HSVColormapFunction<float, RGBA> HSVFunction;

    HSVFunction::Pointer hsvFunc = HSVFunction::New();
    hsvFunc->SetMinimumInputValue(_md[0]);
    hsvFunc->SetMaximumInputValue(_md[1]);
    
    painter->setClipRect(option->exposedRect);
    const int n = m_Counts;
    painter->setPen(m_Pen);
    QPointF center;
    QBrush brush = m_Brush;
    for (int i = 0; i < n; i++) {
        center.setX(m_Particles[i].x[0]);
        center.setY(m_Particles[i].x[1]);
        QColor color;
        if (colorMode == Density) {
            RGBA rgba = hsvFunc->operator()(m_Particles[i].density);
            brush.setColor(QColor(rgba[0], rgba[1], rgba[2]));
        }
        painter->setBrush(brush);
        painter->drawEllipse(center, 1, 1);
    }
}