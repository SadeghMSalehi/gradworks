#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QtGui/QMainWindow>
#include "ui_mainwindow.h"

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

private:
    Ui::MainWindowClass ui;
};

#endif // MAINWINDOW_H
