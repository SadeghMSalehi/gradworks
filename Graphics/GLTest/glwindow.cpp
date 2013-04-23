//
//
// OpenGL Test application
//
//
//

#include "glwindow.h"

GLWindow::GLWindow(QWidget* parent): QMainWindow(parent) {
    ui.setupUi(this);
    ui.graphicsView->setScene(&_scene);
    ui.graphicsView->hide();

    connect(&_timer, SIGNAL(timeout()), SLOT(tickTimer()));
    connect(ui.actionStart, SIGNAL(triggered()), SLOT(startTimer()));
    connect(ui.actionStop, SIGNAL(triggered()), SLOT(stopTimer()));

    ui.widget->setFocus();
}

GLWindow::~GLWindow() {

}

void GLWindow::startTimer() {
    _timer.start(10);
}

void GLWindow::tickTimer() {
    ui.widget->update();
}

void GLWindow::stopTimer() {
    _timer.stop();
}