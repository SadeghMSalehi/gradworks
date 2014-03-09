#include <QApplication>

#include "glwindow.h"


int main(int argc, char* argv[]) {
    QApplication apps(argc, argv);
    GLWindow mainWindow;
    mainWindow.show();
    return apps.exec();
}

