//
//  myQGraphicsView.h
//  ParticlesGUI
//
//  Created by Joohwi Lee on 11/19/12.
//
//

#ifndef __ParticlesGUI__myQGraphicsView__
#define __ParticlesGUI__myQGraphicsView__

#include <iostream>
#include <QGraphicsView>

class myQGraphicsView : public QGraphicsView {
    Q_OBJECT
public:
    myQGraphicsView(QWidget* widget = 0) : QGraphicsView(widget) {
        
    }
    
    virtual ~myQGraphicsView() {
        
    }
    
signals:
    void mouseReleased(QMouseEvent* event);
    void mousePressed(QMouseEvent* event);
    
protected:
    virtual void mouseReleaseEvent(QMouseEvent* event);
    virtual void mousePressEvent(QMouseEvent* event);
};

#endif /* defined(__ParticlesGUI__myQGraphicsView__) */
