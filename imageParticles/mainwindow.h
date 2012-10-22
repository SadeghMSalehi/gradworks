#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtGui/QMainWindow>
#include "QGraphicsScene"
#include "QListWidget"

#include "ui_particlesUI.h"

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = 0);
    ~MainWindow();
    
public slots:
    void on_actionOpen_triggered();
    void on_actionClose_triggered();
    void on_actionZoomIn_triggered();
    void on_actionZoomOut_triggered();
    void on_actionDeploy_triggered();
    void on_actionRun_triggered();
    void on_listWidget_itemClicked(QListWidgetItem * item);
    void on_listWidget_currentRowChanged(int currentRow);

public:
    void updateScene();
    void addImage(QString filename);
    void resizeEvent(QResizeEvent* event);

private:
    Ui::MainWindow ui;
    QGraphicsScene gs;
};

#endif