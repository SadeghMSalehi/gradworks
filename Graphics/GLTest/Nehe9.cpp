//
//  Nehe9.cpp
//  gltest
//
//  Created by Joohwi Lee on 4/23/13.
//
//

#include "Nehe9.h"
#include <GLUT/GLUT.h>
#include <QFileInfo>
#include <iostream>

using namespace std;

void Nehe::initBlending() {
    glShadeModel(GL_SMOOTH);                // Enable Smooth Shading
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);           // Set The Blending Function For Translucency
    glEnable(GL_BLEND);                 // Enable Blending
}

void Nehe::loadTextures() {
    glEnable(GL_TEXTURE_2D);
    glGenTextures(1, _texture);

    for (int i = 0; i < 1; i++) {
        QImage image;
        image.load("/Users/joohwi/Downloads/Star.bmp");
        QImage glImage = QGLWidget::convertToGLFormat(image);

        glBindTexture(GL_TEXTURE_2D, _texture[0]);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, glImage.width(), glImage.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, glImage.scanLine(0));
    }
}

void Nehe::resizeGL(int w, int h) {
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    //    gluOrtho2D(0, w, 0, h); // set origin to bottom left corner
    gluPerspective(45., GLfloat(w)/GLfloat(h), 0.1f, 100.0f);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);  // Really Nice Perspective Calculations
}



void Nehe9::initializeGL() {
    glClearColor(0.0f, 0.0f, 0.0f, 0.5f);           // Black Background
    glClearDepth(1.0f);                 // Depth Buffer Setup

    loadTextures();
    initBlending();

    for (int i = 0; i < 50; i++) {
        stars[i].angle = 0;
        stars[i].dist = (i / 50.0f) * 5.0f;
        stars[i].r = rand() % 256;
        stars[i].g = rand() % 256;
        stars[i].b = rand() % 256;
    }
}


void Nehe9::paintGL() {
    glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);


    glClearColor(0,0,0,1);          // We'll Clear To The Color Of The Fog ( Modified )


    glBindTexture(GL_TEXTURE_2D, _texture[0]);
    for (int i = 0; i < 50; i++) {
        glLoadIdentity();
        glTranslatef(0, 0, zoom);
        glRotatef(tilt,  1.0f, 0, 0);
        glRotatef(stars[i].angle, 0, 1, 0);
        glTranslatef(stars[i].dist, 0, 0);
        glRotatef(-stars[i].angle, 0, 1, 0);
        glRotatef(-tilt, 1, 0, 0);

        glColor4ub(stars[(50-i)-1].r, stars[(50-i)-1].g, stars[(50-i)-1].b, 255);
        glBegin(GL_QUADS);
        glTexCoord2f(0, 0); glVertex2f(-1, -1);
        glTexCoord2f(1, 0); glVertex2f(1, -1);
        glTexCoord2f(1, 1); glVertex2f(1, 1);
        glTexCoord2f(0, 1); glVertex2f(-1, 1);
        glEnd();

        glRotatef(spin, 0, 0, 1);
        glColor4ub(stars[i].r, stars[i].g, stars[i].b, 255);

        glBegin(GL_QUADS);
        glTexCoord2f(0, 0); glVertex2f(-1, -1);
        glTexCoord2f(1, 0); glVertex2f(1, -1);
        glTexCoord2f(1, 1); glVertex2f(1, 1);
        glTexCoord2f(0, 1); glVertex2f(-1, 1);
        glEnd();

        spin += 0.01f;
        stars[i].angle += (float(i) / 50.f);
        stars[i].dist -= 0.01f;

        if (stars[i].dist < 0.f) {
            stars[i].dist += 5.0f;
            stars[i].r = rand() % 256;
            stars[i].g = rand() % 256;
            stars[i].b = rand() % 256;
        }
    }

    _widget->swapBuffers();
}

void Nehe9::keyPressed(int key) {
    switch (key) {
        case Qt::Key_Up:
            tilt -= 0.5f;
            break;
        case Qt::Key_Down:
            tilt += 0.5f;
            break;
        case Qt::Key_I:
            zoom -= 0.2f;
            break;
        case Qt::Key_K:
            zoom += 0.2f;
            break;
    }
    update();
}



CubeGeometry::CubeGeometry() {

}

CubeGeometry::~CubeGeometry() {

}

void CubeGeometry::init() {
    initializeGLFunctions();
    glGenBuffers(2, _vboIds);

    initGeometry();
}

void CubeGeometry::initGeometry() {
    // For cube we would need only 8 vertices but we have to
    // duplicate vertex for each face because texture coordinate
    // is different.
    VertexData vertices[] = {
        // Vertex data for face 0
        {QVector3D(-1.0, -1.0,  1.0), QVector2D(0.0, 0.0)},  // v0
        {QVector3D( 1.0, -1.0,  1.0), QVector2D(0.33, 0.0)}, // v1
        {QVector3D(-1.0,  1.0,  1.0), QVector2D(0.0, 0.5)},  // v2
        {QVector3D( 1.0,  1.0,  1.0), QVector2D(0.33, 0.5)}, // v3

        // Vertex data for face 1
        {QVector3D( 1.0, -1.0,  1.0), QVector2D( 0.0, 0.5)}, // v4
        {QVector3D( 1.0, -1.0, -1.0), QVector2D(0.33, 0.5)}, // v5
        {QVector3D( 1.0,  1.0,  1.0), QVector2D(0.0, 1.0)},  // v6
        {QVector3D( 1.0,  1.0, -1.0), QVector2D(0.33, 1.0)}, // v7

        // Vertex data for face 2
        {QVector3D( 1.0, -1.0, -1.0), QVector2D(0.66, 0.5)}, // v8
        {QVector3D(-1.0, -1.0, -1.0), QVector2D(1.0, 0.5)},  // v9
        {QVector3D( 1.0,  1.0, -1.0), QVector2D(0.66, 1.0)}, // v10
        {QVector3D(-1.0,  1.0, -1.0), QVector2D(1.0, 1.0)},  // v11

        // Vertex data for face 3
        {QVector3D(-1.0, -1.0, -1.0), QVector2D(0.66, 0.0)}, // v12
        {QVector3D(-1.0, -1.0,  1.0), QVector2D(1.0, 0.0)},  // v13
        {QVector3D(-1.0,  1.0, -1.0), QVector2D(0.66, 0.5)}, // v14
        {QVector3D(-1.0,  1.0,  1.0), QVector2D(1.0, 0.5)},  // v15

        // Vertex data for face 4
        {QVector3D(-1.0, -1.0, -1.0), QVector2D(0.33, 0.0)}, // v16
        {QVector3D( 1.0, -1.0, -1.0), QVector2D(0.66, 0.0)}, // v17
        {QVector3D(-1.0, -1.0,  1.0), QVector2D(0.33, 0.5)}, // v18
        {QVector3D( 1.0, -1.0,  1.0), QVector2D(0.66, 0.5)}, // v19

        // Vertex data for face 5
        {QVector3D(-1.0,  1.0,  1.0), QVector2D(0.33, 0.5)}, // v20
        {QVector3D( 1.0,  1.0,  1.0), QVector2D(0.66, 0.5)}, // v21
        {QVector3D(-1.0,  1.0, -1.0), QVector2D(0.33, 1.0)}, // v22
        {QVector3D( 1.0,  1.0, -1.0), QVector2D(0.66, 1.0)}  // v23
    };

    // Indices for drawing cube faces using triangle strips.
    // Triangle strips can be connected by duplicating indices
    // between the strips. If connecting strips have opposite
    // vertex order then last index of the first strip and first
    // index of the second strip needs to be duplicated. If
    // connecting strips have same vertex order then only last
    // index of the first strip needs to be duplicated.
    GLushort indices[] = {
        0,  1,  2,  3,  3,     // Face 0 - triangle strip ( v0,  v1,  v2,  v3)
        4,  4,  5,  6,  7,  7, // Face 1 - triangle strip ( v4,  v5,  v6,  v7)
        8,  8,  9, 10, 11, 11, // Face 2 - triangle strip ( v8,  v9, v10, v11)
        12, 12, 13, 14, 15, 15, // Face 3 - triangle strip (v12, v13, v14, v15)
        16, 16, 17, 18, 19, 19, // Face 4 - triangle strip (v16, v17, v18, v19)
        20, 20, 21, 22, 23      // Face 5 - triangle strip (v20, v21, v22, v23)
    };

    //! [1]
    // Transfer vertex data to VBO 0
    glBindBuffer(GL_ARRAY_BUFFER, _vboIds[0]);
    glBufferData(GL_ARRAY_BUFFER, 24 * sizeof(VertexData), vertices, GL_STATIC_DRAW);

    // Transfer index data to VBO 1
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _vboIds[1]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 34 * sizeof(GLushort), indices, GL_STATIC_DRAW);
    //! [1]

    cout << "Buffer Bind OK: " << _vboIds[0] << ", " << _vboIds[1] << endl;
}

void CubeGeometry::draw(QGLShaderProgram *program) {

    // Tell OpenGL which VBOs to use
    glBindBuffer(GL_ARRAY_BUFFER, _vboIds[0]);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _vboIds[1]);

    // Offset for position
    int offset = 0;

    // Tell OpenGL programmable pipeline how to locate vertex position data
    int vertexLocation = program->attributeLocation("a_position");
    program->enableAttributeArray(vertexLocation);
    glVertexAttribPointer(vertexLocation, 3, GL_FLOAT, GL_FALSE, sizeof(VertexData), (const void *)offset);

    // Offset for texture coordinate
    offset += sizeof(QVector3D);

    // Tell OpenGL programmable pipeline how to locate vertex texture coordinate data
    int texcoordLocation = program->attributeLocation("a_texcoord");
    program->enableAttributeArray(texcoordLocation);
    glVertexAttribPointer(texcoordLocation, 2, GL_FLOAT, GL_FALSE, sizeof(VertexData), (const void *)offset);

    // Draw cube geometry using indices from VBO 1
    glDrawElements(GL_TRIANGLE_STRIP, 34, GL_UNSIGNED_SHORT, 0);

    cout << "GL Drawing ..." << endl;
}

ShaderModule::ShaderModule() {
    _shaderProgram = NULL;
    _vertexShader = NULL;
    _fragmentShader = NULL;
}

ShaderModule::~ShaderModule() {

}


void ShaderModule::initializeGL() {
    //    glClearColor(0.25f, 0.25f, 0.4f, 0.0f);
    glClearColor(0,0,0,1);

    _vshaderFile = "/tools/gradworks/Graphics/GLTest/basic.vsh";
    _fshaderFile = "/tools/gradworks/Graphics/GLTest/basic.fsh";

    loadShader();
    loadTexture();

    // Enable depth buffer
    glEnable(GL_DEPTH_TEST);

    // Enable back face culling
    glEnable(GL_CULL_FACE);

    _geom.init();
}

void ShaderModule::loadTexture() {
    // Loading cube.png to texture unit 0
    glActiveTexture(GL_TEXTURE0);
    glEnable(GL_TEXTURE_2D);

    QImage image;
    image.load("/Users/joohwi/Downloads/Clown-Fish-05.jpg");
    _texture[0] = _widget->bindTexture(image);

    // Set nearest filtering mode for texture minification
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    // Set bilinear filtering mode for texture magnification
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    // Wrap texture coordinates by repeating
    // f.ex. texture coordinate (1.1, 1.2) is same as (0.1, 0.2)
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
}

bool ShaderModule::compileShaderFile(QGLShader::ShaderTypeBit shaderType, QGLShader* &shader, QString& file) {
    if (shader != NULL) {
        delete shader;
        shader = NULL;
    }
    QFileInfo fileInfo(file);
    if (fileInfo.exists()) {
        shader = new QGLShader(shaderType);
        if (shader->compileSourceFile(file)) {
            _shaderProgram->addShader(shader);
        } else {
            qWarning() << file << " compile error!";
            return false;
        }
    } else {
        qWarning() << file << " not found";
        return false;
    }
    return true;
}

void ShaderModule::loadShader() {
    if (_shaderProgram) {
        _shaderProgram->release();
        _shaderProgram->removeAllShaders();
    } else {
        _shaderProgram = new QGLShaderProgram();
    }

    if (compileShaderFile(QGLShader::Vertex, _vertexShader, _vshaderFile)) {
        cout << _vshaderFile.toStdString() << " compiled ok." << endl;
    }

    if (compileShaderFile(QGLShader::Fragment, _fragmentShader, _fshaderFile)) {
        cout << _fshaderFile.toStdString() << " compiled ok." << endl;
    }

    if (!_shaderProgram->link()) {
        qWarning() << "ShaderLinkerError: " << _shaderProgram->log();
    } else {
        _shaderProgram->bind();
        cout << "Program bind ok." << endl;
    }
}



void ShaderModule::resizeGL(int width, int height) {
    glViewport(0, 0, width, height);
    projection.setToIdentity();
    projection.perspective(45., ((GLfloat)width) / ((GLfloat)height), 3.0f, 7.0f);

    qDebug() << projection;
}

void ShaderModule::paintGL() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    QMatrix4x4 matrix;
    matrix.translate(0,0,-5);
    matrix.rotate(rotation);

//    qDebug() << matrix << endl;

    _shaderProgram->setUniformValue("mvp_matrix", matrix);
    _shaderProgram->setUniformValue("texture", 0);

    _geom.draw(_shaderProgram);
}