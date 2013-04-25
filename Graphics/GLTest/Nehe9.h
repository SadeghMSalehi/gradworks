//
//  Nehe9.h
//  gltest
//
//  Created by Joohwi Lee on 4/23/13.
//
//

#ifndef __gltest__Nehe9__
#define __gltest__Nehe9__

#include <iostream>

#include "QGLProxyWidget.h"
#include <QtOpenGL/QGLShaderProgram>
#include <QtOpenGL/QGLShader>
#include <QtOpenGL/QGLFunctions>

class Nehe: public QGLModule {
public:
    virtual void resizeGL(int w, int h);
    virtual void loadTextures();
    virtual void initBlending();

protected:
    GLuint _texture[1];
};

class Nehe9: public Nehe {
private:
    class Star {
    public:
        int r, g, b;
        float dist;
        float angle;
    };

    Star stars[50];

    float zoom;
    float tilt;
    float spin;

public:
    Nehe9(): Nehe() {
        zoom = -15.0f;
        tilt = 90.0f;
        spin = 0.0f;
    }

    virtual void initializeGL();
    virtual void paintGL();
    virtual void keyPressed(int key);
};


class CubeGeometry: protected QGLFunctions {
public:
    CubeGeometry();
    virtual ~CubeGeometry();

    void init();
    void draw(QGLShaderProgram* program);

private:
    void initGeometry();
    
    GLuint _vboIds[2];

    struct VertexData
    {
        QVector3D position;
        QVector2D texCoord;
    };
};

class ShaderModule: public QGLModule {
public:
    ShaderModule();
    virtual ~ShaderModule();
    
    virtual void loadTexture();
    virtual void loadShader();
    virtual bool compileShaderFile(QGLShader::ShaderTypeBit shaderType, QGLShader* &shader, QString& file);

    virtual void initializeGL();
    virtual void resizeGL(int w, int h);
    virtual void paintGL();
private:
    QMatrix4x4 projection;
    QVector2D mousePressPosition;
    QVector3D rotationAxis;
    qreal angularSpeed;
    QQuaternion rotation;

    QGLShader* _vertexShader;
    QGLShader* _fragmentShader;
    QGLShaderProgram* _shaderProgram;
    QString _vshaderFile;
    QString _fshaderFile;

    GLuint _texture[1];
    CubeGeometry _geom;
};

#endif /* defined(__gltest__Nehe9__) */
