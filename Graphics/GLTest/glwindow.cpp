//
//
// OpenGL Test application
//
//
//

#include "glwindow.h"
#include "Nehe9.h"

GLWindow::GLWindow(QWidget* parent): QMainWindow(parent) {
    ui.setupUi(this);
    ui.graphicsView->setScene(&_scene);
    ui.graphicsView->hide();

    connect(&_timer, SIGNAL(timeout()), SLOT(tickTimer()));
    connect(ui.actionStart, SIGNAL(triggered()), SLOT(startTimer()));
    connect(ui.actionStop, SIGNAL(triggered()), SLOT(stopTimer()));

    ui.glWidget->setFocus();

    static ShaderModule module;
    ui.glWidget->setModule(&module);
}

GLWindow::~GLWindow() {

}

void GLWindow::startTimer() {
    _timer.start(10);
}

void GLWindow::tickTimer() {
    ui.glWidget->update();
}

void GLWindow::stopTimer() {
    _timer.stop();
}