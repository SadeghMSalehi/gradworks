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

    void newMenu(void);
    void quitMenu(void);

private:
    Ui::mainwindowClass ui;
};

#endif // MAINWINDOW_H
