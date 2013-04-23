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

struct TexQuad {
    float xcoords[16];
    float tcoords[8];
    
    TexQuad() {
        memset(xcoords, 0, sizeof(xcoords));
        memset(tcoords, 0, sizeof(tcoords));
    }
};

typedef std::vector<TexQuad> TexQuadArray;

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
    
    void loadTextures();
    
private:
    GLuint _texture[1];
};

#endif /* defined(__gltest__QGlyphWidget__) */
