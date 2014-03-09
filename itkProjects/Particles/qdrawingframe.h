//
//  qdrawingframe.h
//  ParticleGuidedRegistration
//
//  Created by Joohwi Lee on 4/3/13.
//
//

#ifndef __ParticleGuidedRegistration__qdrawingframe__
#define __ParticleGuidedRegistration__qdrawingframe__

#include <iostream>
#include <QFrame>
#include <QPainterPath>

class QDrawingFrame: public QFrame {
Q_OBJECT
public:
    QDrawingFrame(QWidget* parent, Qt::WindowType type = Qt::Widget);
    virtual ~QDrawingFrame();

protected:
    virtual void paintEvent(QPaintEvent* event);
    virtual void mousePressEvent(QMouseEvent* event);
    virtual void mouseMoveEvent(QMouseEvent* event);
    virtual void mouseReleaseEvent(QMouseEvent* event);

    QPainterPath _drawingPath;
};

#endif /* defined(__ParticleGuidedRegistration__qdrawingframe__) */
