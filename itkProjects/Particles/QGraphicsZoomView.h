//
//  QGraphicsZoomView.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 5/4/13.
//
//

#ifndef __ParticleGuidedRegistration__QGraphicsZoomView__
#define __ParticleGuidedRegistration__QGraphicsZoomView__

#include <iostream>
#include <QGraphicsView>
#include <QWidget>

class QGraphicsZoomView: public QGraphicsView {
    Q_OBJECT
public:
    QGraphicsZoomView(QWidget* parent = NULL);
    virtual ~QGraphicsZoomView();

    virtual void enterEvent(QEvent* event);
    virtual void leaveEvent(QEvent* event);

    virtual void mouseReleaseEvent(QMouseEvent* event);
    virtual void mouseMoveEvent(QMouseEvent* event);
    virtual void mousePressEvent(QMouseEvent* event);

private:
    QPointF _trackingStartPoint;
    float _scale;
};

#endif /* defined(__ParticleGuidedRegistration__QGraphicsZoomView__) */
