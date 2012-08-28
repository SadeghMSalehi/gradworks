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

private:
    Ui::MainWindowClass ui;
};

#endif // MAINWINDOW_H
