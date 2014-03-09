//
//  QGLProxyWidget.h
//  gltest
//
//  Created by Joohwi Lee on 4/23/13.
//
//

#ifndef __gltest__QGLProxyWidget__
#define __gltest__QGLProxyWidget__

#include <iostream>
#include <QGLWidget>
#include <QKeyEvent>

class QGLModule {
public:
    QGLModule(): _widget(NULL) {}
    virtual ~QGLModule() {}

    inline void setGLWidget(QGLWidget* widget) {
        _widget = widget;
    }

    virtual void initializeGL() = 0;
    virtual void resizeGL(int w, int h) = 0;
    virtual void paintGL() = 0;

    virtual void keyPressed(int key) {}
    virtual void keyReleased(int key) {}
    virtual void mousePressed(QMouseEvent* event) {}
    virtual void mouseMoved(QMouseEvent* event) {}
    virtual void mouseReleased(QMouseEvent* event) {}

    virtual void update() {
        if (_widget) {
            _widget->update();
        }
    }

protected:
    QGLWidget* _widget;
};

class QGLProxyWidget: public QGLWidget {
    Q_OBJECT
    
public:
    QGLProxyWidget(QWidget * parent = 0, const QGLWidget * shareWidget = 0, Qt::WindowFlags f = 0): QGLWidget(parent, shareWidget, f) {
        _glModule = NULL;
    }

    virtual ~QGLProxyWidget() {

    }

    inline void setModule(QGLModule* module) {
        _glModule = module;
        _glModule->setGLWidget(this);
    }

protected:
    inline void initializeGL() {
        if (_glModule != NULL) {
            _glModule->initializeGL();
        }
    }

    inline void resizeGL(int w, int h) {
        if (_glModule != NULL) {
            _glModule->resizeGL(w, h);
        }
    }

    inline void paintGL() {
        if (_glModule != NULL) {
            _glModule->paintGL();
        }
    }

    inline void keyPressEvent(QKeyEvent* event) {
        if (_glModule != NULL) {
            _glModule->keyPressed(event->key());
        }
    }
    inline void keyReleaseEvent(QKeyEvent* event) {
        if (_glModule != NULL) {
            _glModule->keyReleased(event->key());
        }
    }

    inline void mousePressEvent(QMouseEvent* event) {
        if (_glModule != NULL) {
            _glModule->mousePressed(event);
        }
    }

    inline void mouseMoveEvent(QMouseEvent* event) {
        if (_glModule != NULL) {
            _glModule->mouseMoved(event);
        }
    }

    inline void mouseReleaseEvent(QMouseEvent* event) {
        if (_glModule != NULL) {
            _glModule->mouseReleased(event);
        }
    }

private:
    QGLModule* _glModule;
};


#endif /* defined(__gltest__QGLProxyWidget__) */
