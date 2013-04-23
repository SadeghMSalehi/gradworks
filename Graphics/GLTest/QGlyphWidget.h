//
//  QGlyphWidget.h
//  gltest
//
//  Created by Joohwi Lee on 4/22/13.
//
//

#ifndef __gltest__QGlyphWidget__
#define __gltest__QGlyphWidget__

#include <iostream>
#include <QGLWidget>
#include <QKeyEvent>

class QGlyphWidget: public QGLWidget {
    Q_OBJECT
public:
    QGlyphWidget( QWidget * parent = 0, const QGLWidget * shareWidget = 0, Qt::WindowFlags f = 0 );
    virtual ~QGlyphWidget();

protected:
    void initializeGL();
    void resizeGL(int w, int h);
    void paintGL();

    void keyReleaseEvent(QKeyEvent*);
};

#endif /* defined(__gltest__QGlyphWidget__) */
