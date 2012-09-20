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
    void on_action_Open_triggered(bool checked);
    void exit(void);

public:
    void drawImage();

private:
    Ui::MainWindowClass ui;
    QGraphicsScene _scene;
    itkMyCore _core;
};

#endif // MAINWINDOW_H
