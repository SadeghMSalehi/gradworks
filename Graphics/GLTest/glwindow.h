#ifndef GLWINDOW_H
#define GLWINDOW_H

#include <QMainWindow>
#include <QGraphicsScene>
#include <QTimer>

#include "ui_glwindow.h"

class GLWindow: public QMainWindow {
    Q_OBJECT
public:
    GLWindow(QWidget* parent = NULL);
    virtual ~GLWindow();

public slots:
    void startTimer();
    void tickTimer();
    void stopTimer();

private:
    QTimer _timer;
    Ui_GLWindow ui;
    QGraphicsScene _scene;
};

#endif