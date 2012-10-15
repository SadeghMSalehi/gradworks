#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtGui/QMainWindow>
#include "ui_mainwindow.h"
#include <vtkRenderer.h>
#include "PNSMath.h"
#include "MatrixCode.h"

class vtkPolyData;
class vtkGenericOpenGLRenderWindow;
class QVTKInteractor;

class MainWindow: public QMainWindow {
    Q_OBJECT

public:
    MainWindow(QWidget* parent = NULL);
    ~MainWindow();

public slots:
    void on_resetButton_clicked();
    void on_addButton_clicked();
    void on_runButton_clicked();
    void on_testButton_clicked();

private:
    void AddPolyData(vtkPolyData* poly, float r = 1, float g = 0, float b = 0, float opacity = 1.0);
    Ui::MainWindow ui;
    vtkRenderer* m_Renderer;
    QVTKInteractor* m_Interactor;
    PNSMath m_Math;
};

#endif