#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtGui/QMainWindow>
#include "QGraphicsScene"
#include "QListWidget"
#include "QTimer"

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
    void on_actionContinueOptimization_triggered();
    void on_actionPointSave_triggered();
    void on_actionPointLoad_triggered();
    void on_actionPlayTrace_triggered();
    void on_actionSaveTrace_triggered();
    void on_actionLoadTrace_triggered();
    void on_actionShowPlotWindow_triggered();
    void on_actionNBodySimulation_triggered();
    void on_listWidget_itemClicked(QListWidgetItem * item);
    void on_listWidget_currentRowChanged(int currentRow);
    void particleAnimationTimeout();
    void on_graphicsView_mousePressed(QMouseEvent* event);
    void on_animationInterval_valueChanged(int value);
    void on_optiCG_toggled(bool toggled);
    void on_optiLBFG_toggled(bool toggled);
    void on_optiGD_toggled(bool toggled);
    void updateScene();

public:
    void playScene();
    void addImage(QString filename);
    void resizeEvent(QResizeEvent* event);

    void runOptimization();

private:
    Ui::MainWindow ui;
    QGraphicsScene gs;
    QTimer m_Timer;
    QActionGroup* m_particleColors;
};

#endif
