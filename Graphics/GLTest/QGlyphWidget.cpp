//
//  QGlyphWidget.cpp
//  gltest
//
//  Created by Joohwi Lee on 4/22/13.
//
//

#include "QGlyphWidget.h"

#include <GLUT/GLUT.h>

using namespace std;

QGlyphWidget::QGlyphWidget(QWidget* parent, const QGLWidget* shareWidget, Qt::WindowFlags f): QGLWidget(QGLFormat(QGL::DoubleBuffer|QGL::DepthBuffer|QGL::Rgba|QGL::SampleBuffers|QGL::AlphaChannel), parent, shareWidget, f) {
}

QGlyphWidget::~QGlyphWidget() {
}

void QGlyphWidget::initializeGL() {
    glClearColor(0.2f, 0.2f, 0.2f, 0);
    glClearDepth(1.0f);

    glEnable(GL_DEPTH_TEST);
//    glDisable(GL_DEPTH_TEST);
//    glEnable(GL_BLEND);
//    glBlendFunc(GL_SRC_ALPHA,GL_ONE);

    glEnable(GL_TEXTURE_2D);                        // Enable Texture Mapping ( NEW )
    glDisable(GL_CULL_FACE);
    glDepthFunc(GL_LEQUAL);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

    glEnable(GL_POLYGON_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);

//        glEnable(GL_BLEND);
    //    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glShadeModel(GL_SMOOTH);

    GLfloat diffusionLight[4] = { 1, 1, 1, 1 };
    GLfloat ambientLight[4] = { 0.5, 0.5, 0.5, 1};
    GLfloat positionLight[4] = { 0, 0, 2, 1};

    glLightfv(GL_LIGHT1, GL_AMBIENT, ambientLight);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, diffusionLight);
    glLightfv(GL_LIGHT1, GL_POSITION, positionLight);
    glEnable(GL_LIGHT1);

}

void QGlyphWidget::resizeGL(int w, int h) {
    loadTextures();
    
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
//    gluOrtho2D(0, w, 0, h); // set origin to bottom left corner
    gluPerspective(45., GLfloat(w)/GLfloat(h), 0.1f, 100.0f);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void QGlyphWidget::loadTextures() {
    QImage image;
    image.load("/Users/joohwi/Downloads/Star.bmp");
    QImage glImage = convertToGLFormat(image);
    
    glGenTextures(1, _texture);
    for (int i = 0; i < 1; i++) {
        glBindTexture(GL_TEXTURE_2D, _texture[i]);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, glImage.width(), glImage.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, glImage.scanLine(0));
    }
}

void QGlyphWidget::paintGL() {
    static GLfloat rx = 0.0f, ry = 0.0f, rz = 0.0f;
 
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    float cubecoord[6][12] = {
        {-1,-1,1, 1,-1,1, 1,1,1, -1,1,1},
        {-1,-1,-1, 1,-1,-1, 1,1,-1, -1,1,-1},
        {-1,1,-1, -1,1,1, 1,1,1, 1,1,-1},
        {-1,-1,-1, 1,-1,-1, 1,-1,1, -1,-1,1},
        {1,-1,-1, 1,1,-1, 1,1,1, 1,-1,1},
        {-1,-1,-1, -1,-1,1, -1,1,1, -1,1,-1},
    };
    
    float cubetexcoord[6][8] = {
        { 0,0, 1,0, 1,1, 0,1 },
        { 1,0, 1,1, 0,1, 0,0 },
        { 0,1, 0,0, 1,0, 1,1 },
        { 1,1, 0,1, 0,0, 1,0 },
        { 1,1, 0,1, 0,0, 1,0 },
        { 1,0, 1,1, 0,1, 0,0 }
    };
    
    glLoadIdentity();
    glTranslatef(0.0f,0.0f,-5.0f);                      // Move Into The Screen 5 Units

    glRotatef(rx,1,0,0);
    glRotatef(ry,0,1,0);
    glRotatef(rz,0,0,1);
    
    glBindTexture(GL_TEXTURE_2D, _texture[0]);
    glColor4f(1.0f,1.0f,1.0f,0.5f);
    glBegin(GL_QUADS);
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 4; j++) {
            glTexCoord2fv(&cubetexcoord[i][j*2]);
            glVertex3fv(&cubecoord[i][j*3]);
        }
    }
    glEnd();

    rx += 0.2f;
    ry += 0.3f;
    rz += 0.4f;
}

void QGlyphWidget::keyReleaseEvent(QKeyEvent *key) {
    int value = key->key();
    switch (value) {
        case Qt::Key_L:
            glEnable(GL_LIGHTING);
            cout << "Lighting On" << endl;
            break;
        case Qt::Key_K:
            glDisable(GL_LIGHTING);
            cout << "Lighting Off" << endl;
            break;
    }
    update();
}

void QGlyphWidget::mousePressEvent(QMouseEvent *event) {
    event->accept();
    _startingPoint = event->pos();
    _startingPoint.setX(_startingPoint.x()/((width()-1.)/2.) - 1);
    _startingPoint.setY(_startingPoint.y()/((height()-1.)/2.) - 1);
}

void QGlyphWidget::mouseMoveEvent(QMouseEvent *event) {

}

void QGlyphWidget::mouseReleaseEvent(QMouseEvent *event) {

}