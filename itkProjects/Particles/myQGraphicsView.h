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
    enum MouseMode { None, Zooming };

    myQGraphicsView(QWidget* widget = 0) : QGraphicsView(widget) {
        _rightMode = None;
    }
    
    virtual ~myQGraphicsView() {
        
    }
    
signals:
    void mouseReleased(QMouseEvent* event);
    void mousePressed(QMouseEvent* event);
    void mouseMoved(QMouseEvent* event);
    
protected:
    MouseMode _rightMode;
    QPoint _buttonPressPos;
    QRect _buttonPressViewRect;

    virtual void enterEvent(QEvent* event);
    virtual void leaveEvent(QEvent* event);

    virtual void mouseReleaseEvent(QMouseEvent* event);
    virtual void mouseMoveEvent(QMouseEvent* event);
    virtual void mousePressEvent(QMouseEvent* event);
};

#endif /* defined(__ParticlesGUI__myQGraphicsView__) */
