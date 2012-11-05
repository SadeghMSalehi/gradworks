//
//  qmygraphicsview.h
//  imageParticles
//
//  Created by Joohwi Lee on 11/2/12.
//
//

#ifndef __imageParticles__qmygraphicsview__
#define __imageParticles__qmygraphicsview__

#include <iostream>
#include <QGraphicsView>

class QMyGraphicsView : public QGraphicsView {
    Q_OBJECT
public:
    QMyGraphicsView(QWidget* widget = 0) : QGraphicsView(widget) {

    }
    virtual ~QMyGraphicsView() {

    }
signals:
    void mouseReleased(QMouseEvent* event);
    void mousePressed(QMouseEvent* event);
    
protected:
    virtual void mouseReleaseEvent(QMouseEvent* event);
    virtual void mousePressEvent(QMouseEvent* event);
};


#endif /* defined(__imageParticles__qmygraphicsview__) */
