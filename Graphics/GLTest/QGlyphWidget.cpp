//
//  QGlyphWidget.cpp
//  gltest
//
//  Created by Joohwi Lee on 4/22/13.
//
//

#include "QGlyphWidget.h"

#include <GLUT/GLUT.h>

QGlyphWidget::QGlyphWidget(QWidget* parent, const QGLWidget* shareWidget, Qt::WindowFlags f): QGLWidget(QGLFormat(QGL::DoubleBuffer|QGL::DepthBuffer|QGL::Rgba|QGL::SampleBuffers|QGL::AlphaChannel), parent, shareWidget, f) {
}

QGlyphWidget::~QGlyphWidget() {
}

void QGlyphWidget::initializeGL() {
    glClearColor(0.2f, 0.2f, 0.2f, 0);
    glClearDepth(1.0f);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    glEnable(GL_BLEND);
    glEnable(GL_POLYGON_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glShadeModel(GL_SMOOTH);
}

void QGlyphWidget::resizeGL(int w, int h) {
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
//    gluOrtho2D(0, w, 0, h); // set origin to bottom left corner
    gluPerspective(45., GLfloat(w)/GLfloat(h), 0.1f, 100.0f);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void QGlyphWidget::paintGL() {
    static GLfloat rTri = 0.0f;
    static GLfloat rQuad = 0.0f;

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glLoadIdentity();                       // Reset The View
    glTranslatef(-1.5f,0.0f,-6.0f);                 // Move Left 1.5 Units And Into The Screen 6.0
    glRotatef(rTri, 0, 0.f, 1.0f);
    glBegin(GL_TRIANGLES);                      // Drawing Using Triangles
    glColor3f( 1.0f, 0.0f, 0.0f);          // Set The Color To Red
    glVertex3f( 0.0f, 1.0f, 0.0f);              // Top
    glColor3f( 0.0f, 1.0f, 0.0f);          // Set The Color To Green
    glVertex3f(-1.0f,-1.0f, 0.0f);              // Bottom Left
    glColor3f(0.0f,0.0f,1.0f);
    glVertex3f( 1.0f,-1.0f, 0.0f);              // Bottom Right
    glEnd();                            // Finished Drawing The Triangle

    glLoadIdentity();
    glColor3f(1.0f,1.0f,1.0f);
    glTranslatef(3.0f,0.0f,0.0f);                   // Move Right 3 Units
//    glRotatef(rQuad, 0, 1.f, 0);
    glBegin(GL_QUADS);                      // Draw A Quad
    glVertex3f(-1.0f, 1.0f, 0.0f);              // Top Left
    glVertex3f( 1.0f, 1.0f, 0.0f);              // Top Right
    glVertex3f( 1.0f,-1.0f, 0.0f);              // Bottom Right
    glVertex3f(-1.0f,-1.0f, 0.0f);              // Bottom Left
    glEnd();                            // Done Drawing The Quad

    rTri += 0.2f;
}

void QGlyphWidget::keyReleaseEvent(QKeyEvent *key) {
    if (key->key() == ' ') {
        update();
    }
}