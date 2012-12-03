#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtGui/QMainWindow>
#include "ui_mainwindow.h"
#include <vtkRenderer.h>
#include "PNSMath.h"
#include "MatrixCode.h"
#include <map>
#include <string>

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
    void AddPolyData(std::string name, vtkPolyData* poly, float r = 1, float g = 0, float b = 0, float opacity = 1.0);
    void RemovePolyData(std::string name);
    void RemoveAllPolyData();

    Ui::MainWindow ui;
    vtkRenderer* m_Renderer;
    QVTKInteractor* m_Interactor;
    PNSMath m_Math;
    typedef std::map<std::string, vtkProp*> PropMapType;
    typedef std::pair<std::string, vtkProp*> NamedProp;
    PropMapType m_PropMap;

};

#endif
