#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtGui/QMainWindow>
#include "QGraphicsScene"

#include "ui_mainwindow.h"

#include "itkMyCore.h"

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = 0);
    ~MainWindow();

public slots:
    void sayHello(void);
    void on_loadTransformButton_triggered();
    void on_actionOpenSource_triggered();
    void on_actionOpenTarget_triggered();
    void on_actionOpenLabel_triggered();
    void on_actionShowSource_triggered(bool c) {
        if (c) {
            ui.actionShowTarget->setChecked(false);
        }
        drawImage();
    }
    void on_actionShowTarget_triggered(bool c) {
        if (c) {
            ui.actionShowSource->setChecked(false);
        }
        drawImage();
    }
    void on_opacityDial_valueChanged() {
        drawImage();
    }
    void on_actionShowLabel_triggered() {
        drawImage();
    }
    void on_actionZoomIn_triggered() {
        ui.graphicsView->scale(2, 2);
    }
    void on_actionZoomOut_triggered() {
        ui.graphicsView->scale(.5, .5);
    }
    void on_sliceSlider_valueChanged() {
        moveSlice();
    }

    void exit(void);

public:
    void moveSlice();
    void drawImage(bool force = true);

private:
    Ui::MainWindowClass ui;
    QGraphicsScene _scene;
    itkMyCore _core;
    int _currentSlice;
};

#endif // MAINWINDOW_H
